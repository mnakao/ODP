#include "common.h"
static uint64_t *_A_dev, *_B_dev;
static uint64_t *_result, *_result_dev;
static int *_adjacency_dev, *_num_degrees_dev;
static bool _num_degrees_flag = false;
static int _nodes, _degree, _groups, _kind;
static double _mem_usage;

extern "C" void apsp_start_profile();
extern "C" void apsp_end_profile(const char* name, const int kind, const int groups, const double mem_usage, const int procs);
extern "C" double apsp_get_mem_usage(const int kind, const int nodes, const int degree, const int groups,
				     const int *num_degrees, const int procs, const bool is_cpu);
extern "C" int  apsp_get_kind(const int nodes, const int degree, const int* num_degrees, const int groups,
			      const int procs, const bool is_cpu);
extern __global__ void clear_buffers(uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int length);
extern __global__ void popcnt(const uint64_t* __restrict__ B, const int nodes,
                              const unsigned int elements, uint64_t* __restrict__ result);
extern __global__ void matrix_op(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
                                 const int* __restrict__ num_degrees, const int nodes, const int degree, const unsigned int elements, const int based_nodes);
extern __global__ void matrix_op_chunk(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
				       const int* __restrict__ num_degrees, const int nodes, const int degree, const int based_nodes);

static __global__ void init_buffers(uint64_t* __restrict__ A, uint64_t* __restrict__ B,
				    const int nodes, const int groups, const unsigned int elements)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < nodes/groups) {
    unsigned int offset = tid*elements+tid/UINT64_BITS;
    A[offset] = B[offset] = (0x1ULL<<(tid%UINT64_BITS));
    tid += blockDim.x * gridDim.x;
  }
}

static __global__ void init_buffers_saving(uint64_t* __restrict__ A, uint64_t* __restrict__ B,
					   const int nodes, const int groups, const int t)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid<UINT64_BITS*GPU_CHUNK && UINT64_BITS*t*GPU_CHUNK+tid<nodes/groups) {
    unsigned int offset = (UINT64_BITS*t*GPU_CHUNK+tid)*GPU_CHUNK+tid/UINT64_BITS;
    A[offset] = B[offset] = (0x1ULL<<(tid%UINT64_BITS));
    tid += blockDim.x * gridDim.x;
  }
}

static void apsp_cuda_mat(const int* __restrict__ adjacency,
			  int *diameter, long *sum, double *ASPL)
{
  unsigned int elements = (_nodes/_groups+(UINT64_BITS-1))/UINT64_BITS;
  *sum = (long)_nodes * (_nodes - 1);
  *diameter = 1;
  clear_buffers <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes*elements);
  init_buffers  <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes, _groups, elements);

  for(int kk=0;kk<_nodes;kk++){
    matrix_op <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _adjacency_dev, _num_degrees_dev,
				       _nodes, _degree, elements, _nodes/_groups);
    popcnt    <<< BLOCKS, THREADS >>> (_B_dev, _nodes, elements, _result_dev);
    
    cudaMemcpy(_result, _result_dev, sizeof(uint64_t)*BLOCKS, cudaMemcpyDeviceToHost);
    uint64_t num = 0;
    for (int i=0;i<BLOCKS;i++)
      num += _result[i];

     num *= _groups;
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
  unsigned int elements = (_nodes/_groups+UINT64_BITS-1)/UINT64_BITS;
  int parsize = (elements + GPU_CHUNK - 1)/GPU_CHUNK;
  *sum = (long)_nodes * (_nodes - 1);
  *diameter = 1;

  for(int t=0;t<parsize;t++){
    unsigned int kk, l;
    for(l=0; l<UINT64_BITS*GPU_CHUNK && UINT64_BITS*t*GPU_CHUNK+l<_nodes/_groups; l++){}
    clear_buffers       <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes*GPU_CHUNK);
    init_buffers_saving <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes, _groups, t);

    for(kk=0;kk<_nodes;kk++){
      matrix_op_chunk <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _adjacency_dev, _num_degrees_dev,
					       _nodes, _degree, _nodes/_groups);
      popcnt    <<< BLOCKS, THREADS >>> (_B_dev, _nodes, GPU_CHUNK, _result_dev);

      cudaMemcpy(_result, _result_dev, sizeof(uint64_t)*BLOCKS, cudaMemcpyDeviceToHost);
      uint64_t num = 0;
      for (int i=0;i<BLOCKS;i++)
        num += _result[i];

      if(num == (uint64_t)_nodes*l) break;

      // swap A <-> B
      uint64_t* tmp = _A_dev;
      _A_dev = _B_dev;
      _B_dev = tmp;
      
      *sum += ((long)_nodes * l - num) * _groups;
    }
    *diameter = MAX(*diameter, kk+1);
  }

  *ASPL = *sum / (((double)_nodes-1)*_nodes);
  *sum /= 2.0;
}

extern "C" void apsp_cuda_init_s(const int nodes, const int degree,
				 const int* __restrict__ num_degrees, const int groups)
{
  cuInit(0);

  if(nodes % groups != 0)
    ERROR("nodes(%d) must be divisible by group(%d)\n", nodes, groups);

  _kind = apsp_get_kind(nodes, degree, num_degrees, groups, 1, false);
  _mem_usage = apsp_get_mem_usage(_kind, nodes, degree, groups, num_degrees, 1, false);
  size_t s = (_kind == APSP_NORMAL)? (nodes/groups+(UINT64_BITS-1))/UINT64_BITS : GPU_CHUNK;
  s *= nodes * sizeof(uint64_t);

  _nodes = nodes;
  _degree = degree;
  _groups = groups;
  
  cudaMalloc((void**)&_A_dev, s);
  cudaMalloc((void**)&_B_dev, s);
  cudaHostAlloc((void**)&_result,     sizeof(uint64_t)*BLOCKS, cudaHostAllocDefault);
  cudaMalloc((void**)&_result_dev,    sizeof(uint64_t)*BLOCKS);
  cudaMalloc((void**)&_adjacency_dev, sizeof(int)*(nodes/groups)*degree);
  if(num_degrees){
    cudaMalloc((void**)&_num_degrees_dev, sizeof(int)*nodes);
    cudaMemcpy(_num_degrees_dev, num_degrees, sizeof(int)*nodes, cudaMemcpyHostToDevice);
    _num_degrees_flag = true;
  }
}

extern "C" void apsp_cuda_init(const int nodes, const int degree, const int* __restrict__ num_degrees)
{
  apsp_cuda_init_s(nodes, degree, num_degrees, 1);
}

extern "C" void apsp_cuda_finalize()
{
  cudaFree(_A_dev);
  cudaFree(_B_dev);
  cudaFreeHost(_result);
  cudaFree(_result_dev);
  cudaFree(_adjacency_dev);
  if(_num_degrees_flag)
    cudaFree(_num_degrees_dev);
}

extern "C" void apsp_cuda_run(const int* __restrict__ adjacency,
			      int *diameter, long *sum, double *ASPL)
{
  apsp_start_profile();
  cudaMemcpy(_adjacency_dev, adjacency, sizeof(int)*(_nodes/_groups)*_degree, cudaMemcpyHostToDevice);
  
  if(_kind == APSP_NORMAL)
    apsp_cuda_mat       (adjacency, diameter, sum, ASPL);
  else
    apsp_cuda_mat_saving(adjacency, diameter, sum, ASPL);

  if(*diameter > _nodes)
    ERROR("This graph is not connected graph.\n");

  apsp_end_profile("CUDA", _kind, _groups, _mem_usage, 1);
}
