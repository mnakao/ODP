#include "common.h"
static uint64_t *_A_dev, *_B_dev;
static uint64_t *_result, *_result_dev;
static int *_adjacency_dev, *_num_degrees_dev;
static bool _num_degrees_flag = false, _is_profile;
static int _nodes, _degree, _symmetries, _kind;
static double _mem_usage, _elapsed_time;
static unsigned int _times;

extern "C" bool ODP_Check_profile();
extern "C" double ODP_Get_time();
extern "C" void ODP_Profile(const char* name, const int kind, const int symmetries, const double mem_usage,
			    const double elapsed_time, const unsigned int times, const int procs);
extern "C" double ODP_Get_mem_usage(const int kind, const int nodes, const int degree, const int symmetries,
				    const int *num_degrees, const int procs, const bool is_cpu);
extern "C" int  ODP_Get_kind(const int nodes, const int degree, const int* num_degrees, const int symmetries,
			     const int procs, const bool is_cpu);
extern __global__ void ODP_Clear_buffers(uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int length);
extern __global__ void ODP_Popcnt(const uint64_t* __restrict__ B, const int nodes,
				  const unsigned int elements, uint64_t* __restrict__ result);
extern __global__ void ODP_Matmul_cuda(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
				       const int* __restrict__ num_degrees, const int nodes, const int degree, const unsigned int elements, const int based_nodes);
extern __global__ void ODP_Matmul_CHUNK_cuda(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
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

static void aspl_cuda_mat(const int* __restrict__ adjacency,
			  int *diameter, long *sum, double *ASPL)
{
  unsigned int elements = (_nodes/_symmetries+(UINT64_BITS-1))/UINT64_BITS;
  *sum = (long)_nodes * (_nodes - 1);
  *diameter = 1;
  ODP_Clear_buffers <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes*elements);
  init_buffers      <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes, _symmetries, elements);

  for(int kk=0;kk<_nodes;kk++){
    ODP_Matmul_cuda <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _adjacency_dev, _num_degrees_dev,
					      _nodes, _degree, elements, _nodes/_symmetries);
    ODP_Popcnt      <<< BLOCKS, THREADS >>> (_B_dev, _nodes, elements, _result_dev);
    
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

static void aspl_cuda_mat_saving(const int* __restrict__ adjacency,
				 int *diameter, long *sum, double *ASPL)
{
  unsigned int elements = (_nodes/_symmetries+UINT64_BITS-1)/UINT64_BITS;
  int parsize = (elements + GPU_CHUNK - 1)/GPU_CHUNK;
  *sum = (long)_nodes * (_nodes - 1);
  *diameter = 1;

  for(int t=0;t<parsize;t++){
    unsigned int kk, l;
    for(l=0; l<UINT64_BITS*GPU_CHUNK && UINT64_BITS*t*GPU_CHUNK+l<_nodes/_symmetries; l++){}
    ODP_Clear_buffers   <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes*GPU_CHUNK);
    init_buffers_saving <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes, _symmetries, t);

    for(kk=0;kk<_nodes;kk++){
      ODP_Matmul_CHUNK_cuda <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _adjacency_dev, _num_degrees_dev,
						      _nodes, _degree, _nodes/_symmetries);
      ODP_Popcnt            <<< BLOCKS, THREADS >>> (_B_dev, _nodes, GPU_CHUNK, _result_dev);

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

extern "C" void ODP_Init_aspl_cuda_s(const int nodes, const int degree,
				     const int* __restrict__ num_degrees, const int symmetries)
{
  cuInit(0);

  if(nodes % symmetries != 0)
    ERROR("nodes(%d) must be divisible by symmetries(%d)\n", nodes, symmetries);

  _kind = ODP_Get_kind(nodes, degree, num_degrees, symmetries, 1, false);
  _mem_usage = ODP_Get_mem_usage(_kind, nodes, degree, symmetries, num_degrees, 1, false);
  size_t s = (_kind == ASPL_NORMAL)? (nodes/symmetries+(UINT64_BITS-1))/UINT64_BITS : GPU_CHUNK;
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
  _is_profile = ODP_Check_profile();
  _elapsed_time = 0;
  _times = 0;
}

extern "C" void ODP_Init_aspl_cuda(const int nodes, const int degree, const int* __restrict__ num_degrees)
{
  ODP_Init_aspl_cuda_s(nodes, degree, num_degrees, 1);
}

extern "C" void ODP_Finalize_aspl_cuda()
{
  cudaFree(_A_dev);
  cudaFree(_B_dev);
  cudaFreeHost(_result);
  cudaFree(_result_dev);
  cudaFree(_adjacency_dev);
  if(_num_degrees_flag)
    cudaFree(_num_degrees_dev);

  if(_is_profile)
    ODP_Profile("CUDA", _kind, _symmetries, _mem_usage,
		_elapsed_time, _times, 1);
}

extern "C" void ODP_Set_aspl_cuda(const int* __restrict__ adjacency,
			      int *diameter, long *sum, double *ASPL)
{
  double t = ODP_Get_time();
  
  cudaMemcpy(_adjacency_dev, adjacency, sizeof(int)*(_nodes/_symmetries)*_degree, cudaMemcpyHostToDevice);
  
  if(_kind == ASPL_NORMAL)
    aspl_cuda_mat       (adjacency, diameter, sum, ASPL);
  else
    aspl_cuda_mat_saving(adjacency, diameter, sum, ASPL);

  _elapsed_time += ODP_Get_time() - t;
    
  if(*diameter > _nodes)
    ERROR("This graph is not connected graph.\n");

  _times++;
}