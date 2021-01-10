#include "common.h"
static uint64_t *_A_dev, *_B_dev;
static uint64_t *_result, *_result_dev;
static int *_adjacency_dev, *_num_degrees_dev;
static bool _num_degrees_flag = false, _is_profile;
static int _nodes, _degree, _symmetries, _kind;
static double _mem_usage, _elapsed_time;
static unsigned int _times;

extern "C" bool apsp_check_profile();
extern "C" double apsp_get_time();
extern "C" void apsp_profile(const char* name, const int kind, const int symmetries, const double mem_usage,
			     const double elapsed_time, const unsigned int times, const int procs);
extern "C" double apsp_get_mem_usage(const int kind, const int nodes, const int degree, const int symmetries,
				     const int *num_degrees, const int procs, const bool is_cpu);
extern "C" int  apsp_get_kind(const int nodes, const int degree, const int* num_degrees, const int symmetries,
			      const int procs, const bool is_cpu);
extern __global__ void apsp_clear_buffers(uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int length);
extern __global__ void apsp_popcnt(const uint64_t* __restrict__ B, const int nodes,
				   const unsigned int elements, uint64_t* __restrict__ result);
extern __global__ void apsp_matmul_cuda(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
					const int* __restrict__ num_degrees, const int nodes, const int degree, const unsigned int elements, const int based_nodes);
extern __global__ void apsp_matmul_CHUNK_cuda(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
					      const int* __restrict__ num_degrees, const int nodes, const int degree, const int based_nodes);

static __global__ void init_buffers(uint64_t* __restrict__ A, uint64_t* __restrict__ B,
				    const int nodes, const int symmetries, const unsigned int elements)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < nodes/symmetries) {
    unsigned int offset = tid*elements+tid/UINT64_BITS;
    A[offset] = B[offset] = (0x1ULL<<(tid%UINT64_BITS));
    tid += blockDim.x * gridDim.x;
  }
}

static __global__ void init_buffers_saving(uint64_t* __restrict__ A, uint64_t* __restrict__ B,
					   const int nodes, const int symmetries, const int t)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid<UINT64_BITS*GPU_CHUNK && UINT64_BITS*t*GPU_CHUNK+tid<nodes/symmetries) {
    unsigned int offset = (UINT64_BITS*t*GPU_CHUNK+tid)*GPU_CHUNK+tid/UINT64_BITS;
    A[offset] = B[offset] = (0x1ULL<<(tid%UINT64_BITS));
    tid += blockDim.x * gridDim.x;
  }
}

static void apsp_cuda_mat(const int* __restrict__ adjacency,
			  int *diameter, long *sum, double *ASPL)
{
  unsigned int elements = (_nodes/_symmetries+(UINT64_BITS-1))/UINT64_BITS;
  *sum = (long)_nodes * (_nodes - 1);
  *diameter = 1;
  apsp_clear_buffers <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes*elements);
  init_buffers       <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes, _symmetries, elements);

  for(int kk=0;kk<_nodes;kk++){
    apsp_matmul_cuda <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _adjacency_dev, _num_degrees_dev,
					      _nodes, _degree, elements, _nodes/_symmetries);
    apsp_popcnt      <<< BLOCKS, THREADS >>> (_B_dev, _nodes, elements, _result_dev);
    
    cudaMemcpy(_result, _result_dev, sizeof(uint64_t)*BLOCKS, cudaMemcpyDeviceToHost);
    uint64_t num = 0;
    for (int i=0;i<BLOCKS;i++)
      num += _result[i];

    num *= _symmetries;
    if(num == (uint64_t)_nodes*_nodes) break;

    // swap A <-> B
    uint64_t* tmp = _A_dev;
    _A_dev = _B_dev;
    _B_dev = tmp;

    *sum += (long)_nodes * _nodes - num;
    (*diameter) += 1;
  }
  
  *ASPL = *sum / (((double)_nodes-1)*_nodes);
  *sum /= 2.0;
}

static void apsp_cuda_mat_saving(const int* __restrict__ adjacency,
				 int *diameter, long *sum, double *ASPL)
{
  unsigned int elements = (_nodes/_symmetries+UINT64_BITS-1)/UINT64_BITS;
  int parsize = (elements + GPU_CHUNK - 1)/GPU_CHUNK;
  *sum = (long)_nodes * (_nodes - 1);
  *diameter = 1;

  for(int t=0;t<parsize;t++){
    unsigned int kk, l;
    for(l=0; l<UINT64_BITS*GPU_CHUNK && UINT64_BITS*t*GPU_CHUNK+l<_nodes/_symmetries; l++){}
    apsp_clear_buffers  <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes*GPU_CHUNK);
    init_buffers_saving <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes, _symmetries, t);

    for(kk=0;kk<_nodes;kk++){
      apsp_matmul_CHUNK_cuda <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _adjacency_dev, _num_degrees_dev,
						      _nodes, _degree, _nodes/_symmetries);
      apsp_popcnt            <<< BLOCKS, THREADS >>> (_B_dev, _nodes, GPU_CHUNK, _result_dev);

      cudaMemcpy(_result, _result_dev, sizeof(uint64_t)*BLOCKS, cudaMemcpyDeviceToHost);
      uint64_t num = 0;
      for (int i=0;i<BLOCKS;i++)
        num += _result[i];

      if(num == (uint64_t)_nodes*l) break;

      // swap A <-> B
      uint64_t* tmp = _A_dev;
      _A_dev = _B_dev;
      _B_dev = tmp;
      
      *sum += ((long)_nodes * l - num) * _symmetries;
    }
    *diameter = MAX(*diameter, kk+1);
  }

  *ASPL = *sum / (((double)_nodes-1)*_nodes);
  *sum /= 2.0;
}

extern "C" void apsp_cuda_run_init_s(const int nodes, const int degree,
				     const int* __restrict__ num_degrees, const int symmetries)
{
  cuInit(0);

  if(nodes % symmetries != 0)
    ERROR("nodes(%d) must be divisible by symmetries(%d)\n", nodes, symmetries);

  _kind = apsp_get_kind(nodes, degree, num_degrees, symmetries, 1, false);
  _mem_usage = apsp_get_mem_usage(_kind, nodes, degree, symmetries, num_degrees, 1, false);
  size_t s = (_kind == APSP_NORMAL)? (nodes/symmetries+(UINT64_BITS-1))/UINT64_BITS : GPU_CHUNK;
  s *= nodes * sizeof(uint64_t);

  _nodes = nodes;
  _degree = degree;
  _symmetries = symmetries;
  
  cudaMalloc((void**)&_A_dev, s);
  cudaMalloc((void**)&_B_dev, s);
  cudaHostAlloc((void**)&_result,     sizeof(uint64_t)*BLOCKS, cudaHostAllocDefault);
  cudaMalloc((void**)&_result_dev,    sizeof(uint64_t)*BLOCKS);
  cudaMalloc((void**)&_adjacency_dev, sizeof(int)*(nodes/symmetries)*degree);
  if(num_degrees){
    cudaMalloc((void**)&_num_degrees_dev, sizeof(int)*nodes);
    cudaMemcpy(_num_degrees_dev, num_degrees, sizeof(int)*nodes, cudaMemcpyHostToDevice);
    _num_degrees_flag = true;
  }
  _is_profile = apsp_check_profile();
  _elapsed_time = 0;
  _times = 0;
}

extern "C" void apsp_cuda_run_init(const int nodes, const int degree, const int* __restrict__ num_degrees)
{
  apsp_cuda_run_init_s(nodes, degree, num_degrees, 1);
}

extern "C" void apsp_cuda_run_finalize()
{
  cudaFree(_A_dev);
  cudaFree(_B_dev);
  cudaFreeHost(_result);
  cudaFree(_result_dev);
  cudaFree(_adjacency_dev);
  if(_num_degrees_flag)
    cudaFree(_num_degrees_dev);

  if(_is_profile)
    apsp_profile("CUDA", _kind, _symmetries, _mem_usage,
		 _elapsed_time, _times, 1);
}

extern "C" void apsp_cuda_run(const int* __restrict__ adjacency,
			      int *diameter, long *sum, double *ASPL)
{
  double t = apsp_get_time();
  
  cudaMemcpy(_adjacency_dev, adjacency, sizeof(int)*(_nodes/_symmetries)*_degree, cudaMemcpyHostToDevice);
  
  if(_kind == APSP_NORMAL)
    apsp_cuda_mat       (adjacency, diameter, sum, ASPL);
  else
    apsp_cuda_mat_saving(adjacency, diameter, sum, ASPL);

  _elapsed_time += apsp_get_time() - t;
    
  if(*diameter > _nodes)
    ERROR("This graph is not connected graph.\n");

  _times++;
}
