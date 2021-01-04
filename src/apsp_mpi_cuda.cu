#include "common.h"
#include <mpi.h>
static uint64_t *_A_dev, *_B_dev;
static uint64_t *_result, *_result_dev;
static int *_adjacency_dev, *_num_degrees_dev;
static bool _num_degrees_flag = false;
static int _nodes, _degree, _groups, _rank, _procs, _kind;
static double _mem_usage;
static MPI_Comm _comm;

extern __global__ void clear_buffers(uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int length);
extern __global__ void popcnt(const uint64_t* __restrict__ B, const int nodes,
                              const unsigned int elements, uint64_t* __restrict__ result);
extern __global__ void matrix_op(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
                                 const int* __restrict__ num_degrees, const int nodes, const int degree, const unsigned int elements, const int based_nodes);
extern __global__ void matrix_op_chunk(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
                                       const int* __restrict__ num_degrees, const int nodes, const int degree, const int based_nodes);

extern "C" void apsp_start_profile();
extern "C" void apsp_end_profile(const char* name, const int kind, const int groups, const double mem_usage, const int procs);
extern "C" int  apsp_get_kind(const int nodes, const int degree, const int* num_degrees, const int groups,
			      const int procs, const bool is_cpu);
extern "C" double apsp_get_mem_usage(const int kind, const int nodes, const int degree, const int groups,
				     const int *num_degrees, const int procs, const bool is_cpu);

static __global__ void init_buffers(uint64_t* __restrict__ A, uint64_t* __restrict__ B,
				    const int nodes, const int groups, const int t, const int chunk)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid<UINT64_BITS*chunk && UINT64_BITS*t*chunk+tid<nodes/groups) {
    unsigned int offset = (UINT64_BITS*t*chunk+tid)*chunk+tid/UINT64_BITS;
    A[offset] = B[offset] = (0x1ULL<<(tid%UINT64_BITS));
    tid += blockDim.x * gridDim.x;
  }
}

static void apsp_mpi_cuda_mat(const int* __restrict__ adjacency, 
			      int *diameter, long *sum, double *ASPL)
{
  unsigned int elements = (_nodes/_groups+(UINT64_BITS-1))/UINT64_BITS;
  unsigned int chunk = (elements+(_procs-1))/_procs;
  int parsize = (elements+(chunk-1))/chunk;

  *sum = 0.0;
  *diameter = 1;

  for(int t=_rank;t<parsize;t+=_procs){
    unsigned int kk, l;
    for(l=0; l<UINT64_BITS*chunk && UINT64_BITS*t*chunk+l<_nodes/_groups; l++){}
    clear_buffers <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes*chunk);
    init_buffers  <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes, _groups, t, chunk);

    for(kk=0;kk<_nodes;kk++){
      matrix_op <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _adjacency_dev, _num_degrees_dev,
					 _nodes, _degree, chunk, _nodes/_groups);
      popcnt    <<< BLOCKS, THREADS >>> (_B_dev, _nodes, chunk, _result_dev);

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
  
  MPI_Allreduce(MPI_IN_PLACE, diameter, 1, MPI_INT,  MPI_MAX, _comm);
  MPI_Allreduce(MPI_IN_PLACE, sum,      1, MPI_LONG, MPI_SUM, _comm);
  *sum += (long)_nodes * (_nodes - 1);

  *ASPL = *sum / (((double)_nodes-1)*_nodes);
  *sum /= 2.0;
}

static void apsp_mpi_cuda_mat_saving(const int* __restrict__ adjacency,
				     int *diameter, long *sum, double *ASPL)
{
  unsigned int elements = (_nodes/_groups+(UINT64_BITS-1))/UINT64_BITS;
  int parsize = (elements+(GPU_CHUNK-1))/GPU_CHUNK;

  *sum = 0;
  *diameter = 1;

  for(int t=_rank;t<parsize;t+=_procs){
    unsigned int kk, l;
    for(l=0; l<UINT64_BITS*GPU_CHUNK && UINT64_BITS*t*GPU_CHUNK+l<_nodes/_groups; l++){}
    clear_buffers <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes*GPU_CHUNK);
    init_buffers  <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes, _groups, t, GPU_CHUNK);

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

  MPI_Allreduce(MPI_IN_PLACE, diameter, 1, MPI_INT,  MPI_MAX, _comm);
  MPI_Allreduce(MPI_IN_PLACE, sum,      1, MPI_LONG, MPI_SUM, _comm);
  *sum += (long)_nodes * (_nodes - 1);
  
  *ASPL = *sum / (((double)_nodes-1)*_nodes);
  *sum /= 2.0;
}

extern "C" void apsp_mpi_cuda_init_s(const int nodes, const int degree,
				     const int* __restrict__ num_degrees, MPI_Comm comm, const int groups)
{
  cuInit(0);
  
  if(nodes % groups != 0)
    ERROR("nodes(%d) must be divisible by group(%d)\n", nodes, groups);

  MPI_Comm_rank(comm, &_rank);
  MPI_Comm_size(comm, &_procs);
  _kind = apsp_get_kind(nodes, degree, num_degrees, groups, _procs, false);
  _mem_usage = apsp_get_mem_usage(_kind, nodes, degree, groups, num_degrees, _procs, false);
  size_t elements = (nodes/groups+(UINT64_BITS-1))/UINT64_BITS;
  size_t s = (_kind == APSP_NORMAL)? (elements+(_procs-1))/_procs : GPU_CHUNK;
  s *= nodes * sizeof(uint64_t);
  
   _nodes = nodes;
   _degree = degree;
   _groups = groups;
   _comm = comm;
    
  cudaMalloc((void**)&_A_dev, s);
  cudaMalloc((void**)&_B_dev, s);
  cudaHostAlloc((void**)&_result, BLOCKS*sizeof(uint64_t), cudaHostAllocDefault);
  cudaMalloc((void**)&_result_dev,      sizeof(uint64_t)*BLOCKS);
  cudaMalloc((void**)&_adjacency_dev,   sizeof(int)*(nodes/groups)*degree);
  if(num_degrees){
    cudaMalloc((void**)&_num_degrees_dev, sizeof(int)*nodes);
    cudaMemcpy(_num_degrees_dev, num_degrees, sizeof(int)*nodes, cudaMemcpyHostToDevice);
    _num_degrees_flag = true;
  }
}

extern "C" void apsp_mpi_cuda_init(const int nodes, const int degree,
				   const int* __restrict__ num_degrees, MPI_Comm comm)
{
  apsp_mpi_cuda_init_s(nodes, degree, num_degrees, comm, 1);
}

extern "C" void apsp_mpi_cuda_finalize()
{
  cudaFree(_A_dev);
  cudaFree(_B_dev);
  cudaFreeHost(_result);
  cudaFree(_result_dev);
  cudaFree(_adjacency_dev);
  if(_num_degrees_flag)
    cudaFree(_num_degrees_dev);
}

extern "C" void apsp_mpi_cuda_run(const int* __restrict__ adjacency,
				  int *diameter, long *sum, double *ASPL)
{
  MPI_Barrier(_comm);
  if(_rank == 0)
    apsp_start_profile();
  
  cudaMemcpy(_adjacency_dev, adjacency, sizeof(int)*(_nodes/_groups)*_degree, cudaMemcpyHostToDevice);
  
  if(_kind == APSP_NORMAL)
    apsp_mpi_cuda_mat       (adjacency, diameter, sum, ASPL);
  else
    apsp_mpi_cuda_mat_saving(adjacency, diameter, sum, ASPL);

  if(*diameter > _nodes)
    ERROR("This graph is not connected graph.\n");

  if(_rank == 0)
    apsp_end_profile("MPI+CUDA", _kind, _groups, _mem_usage, _procs);
}
