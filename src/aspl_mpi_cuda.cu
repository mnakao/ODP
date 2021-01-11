#include "common.h"
#include <mpi.h>
static uint64_t *_A_dev, *_B_dev;
static uint64_t *_result, *_result_dev;
static int *_adjacency_dev, *_num_degrees_dev;
static bool _num_degrees_flag = false, _is_profile;
static int _nodes, _degree, _symmetries, _rank, _procs, _kind;
static double _mem_usage, _elapsed_time;
static unsigned int _times;
static MPI_Comm _comm;

extern "C" bool ODP_Check_profile();
extern "C" double ODP_Get_time();
extern "C" void ODP_Profile(const char* name, const int kind, const int symmetries, const double mem_usage,
			    const double elapsed_time, const unsigned int times, const int procs);
extern "C" int  ODP_Get_kind(const int nodes, const int degree, const int* num_degrees, const int symmetries,
			     const int procs, const bool is_cpu);
extern "C" double ODP_Get_mem_usage(const int kind, const int nodes, const int degree, const int symmetries,
				    const int *num_degrees, const int procs, const bool is_cpu);
extern __global__ void ODP_Clear_buffers(uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int length);
extern __global__ void ODP_Popcnt(const uint64_t* __restrict__ B, const int nodes,
				  const unsigned int elements, uint64_t* __restrict__ result);
extern __global__ void ODP_Matmul_cuda(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
				       const int* __restrict__ num_degrees, const int nodes, const int degree, const unsigned int elements, const int based_nodes);
extern __global__ void ODP_Matmul_CHUNK_cuda(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
					     const int* __restrict__ num_degrees, const int nodes, const int degree, const int based_nodes);

static __global__ void init_buffers(uint64_t* __restrict__ A, uint64_t* __restrict__ B,
				    const int nodes, const int symmetries, const int t, const int chunk)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid<UINT64_BITS*chunk && UINT64_BITS*t*chunk+tid<nodes/symmetries) {
    unsigned int offset = (UINT64_BITS*t*chunk+tid)*chunk+tid/UINT64_BITS;
    A[offset] = B[offset] = (0x1ULL<<(tid%UINT64_BITS));
    tid += blockDim.x * gridDim.x;
  }
}

static void aspl_mpi_cuda_mat(const int* __restrict__ adjacency, 
			      int *diameter, long *sum, double *ASPL)
{
  unsigned int elements = (_nodes/_symmetries+(UINT64_BITS-1))/UINT64_BITS;
  unsigned int chunk = (elements+(_procs-1))/_procs;
  int parsize = (elements+(chunk-1))/chunk;

  *sum = 0.0;
  *diameter = 1;

  for(int t=_rank;t<parsize;t+=_procs){
    unsigned int kk, l;
    for(l=0; l<UINT64_BITS*chunk && UINT64_BITS*t*chunk+l<_nodes/_symmetries; l++){}
    ODP_Clear_buffers <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes*chunk);
    init_buffers      <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes, _symmetries, t, chunk);

    for(kk=0;kk<_nodes;kk++){
      ODP_Matmul_cuda <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _adjacency_dev, _num_degrees_dev,
						_nodes, _degree, chunk, _nodes/_symmetries);
      ODP_Popcnt      <<< BLOCKS, THREADS >>> (_B_dev, _nodes, chunk, _result_dev);

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
  
  MPI_Allreduce(MPI_IN_PLACE, diameter, 1, MPI_INT,  MPI_MAX, _comm);
  MPI_Allreduce(MPI_IN_PLACE, sum,      1, MPI_LONG, MPI_SUM, _comm);
  *sum += (long)_nodes * (_nodes - 1);

  *ASPL = *sum / (((double)_nodes-1)*_nodes);
  *sum /= 2.0;
}

static void aspl_mpi_cuda_mat_saving(const int* __restrict__ adjacency,
				     int *diameter, long *sum, double *ASPL)
{
  unsigned int elements = (_nodes/_symmetries+(UINT64_BITS-1))/UINT64_BITS;
  int parsize = (elements+(GPU_CHUNK-1))/GPU_CHUNK;

  *sum = 0;
  *diameter = 1;

  for(int t=_rank;t<parsize;t+=_procs){
    unsigned int kk, l;
    for(l=0; l<UINT64_BITS*GPU_CHUNK && UINT64_BITS*t*GPU_CHUNK+l<_nodes/_symmetries; l++){}
    ODP_Clear_buffers <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes*GPU_CHUNK);
    init_buffers      <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes, _symmetries, t, GPU_CHUNK);

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

  MPI_Allreduce(MPI_IN_PLACE, diameter, 1, MPI_INT,  MPI_MAX, _comm);
  MPI_Allreduce(MPI_IN_PLACE, sum,      1, MPI_LONG, MPI_SUM, _comm);
  *sum += (long)_nodes * (_nodes - 1);
  
  *ASPL = *sum / (((double)_nodes-1)*_nodes);
  *sum /= 2.0;
}

extern "C" void ODP_Init_aspl_mpi_cuda_s(const int nodes, const int degree,
					 const int* __restrict__ num_degrees, MPI_Comm comm, const int symmetries)
{
  cuInit(0);
  
  if(nodes % symmetries != 0)
    ERROR("nodes(%d) must be divisible by symmetries(%d)\n", nodes, symmetries);

  MPI_Comm_rank(comm, &_rank);
  MPI_Comm_size(comm, &_procs);
  _kind = ODP_Get_kind(nodes, degree, num_degrees, symmetries, _procs, false);
  _mem_usage = ODP_Get_mem_usage(_kind, nodes, degree, symmetries, num_degrees, _procs, false);
  size_t elements = (nodes/symmetries+(UINT64_BITS-1))/UINT64_BITS;
  size_t s = (_kind == ASPL_NORMAL)? (elements+(_procs-1))/_procs : GPU_CHUNK;
  s *= nodes * sizeof(uint64_t);
  
   _nodes = nodes;
   _degree = degree;
   _symmetries = symmetries;
   _comm = comm;
    
  cudaMalloc((void**)&_A_dev, s);
  cudaMalloc((void**)&_B_dev, s);
  cudaHostAlloc((void**)&_result, BLOCKS*sizeof(uint64_t), cudaHostAllocDefault);
  cudaMalloc((void**)&_result_dev,      sizeof(uint64_t)*BLOCKS);
  cudaMalloc((void**)&_adjacency_dev,   sizeof(int)*(nodes/symmetries)*degree);
  if(num_degrees){
    cudaMalloc((void**)&_num_degrees_dev, sizeof(int)*nodes);
    cudaMemcpy(_num_degrees_dev, num_degrees, sizeof(int)*nodes, cudaMemcpyHostToDevice);
    _num_degrees_flag = true;
  }
  _is_profile = ODP_Check_profile();
  _elapsed_time = 0;
  _times = 0;
}

extern "C" void ODP_Init_aspl_mpi_cuda(const int nodes, const int degree,
				       const int* __restrict__ num_degrees, MPI_Comm comm)
{
  ODP_Init_aspl_mpi_cuda_s(nodes, degree, num_degrees, comm, 1);
}

extern "C" void ODP_Finalize_aspl_mpi_cuda()
{
  cudaFree(_A_dev);
  cudaFree(_B_dev);
  cudaFreeHost(_result);
  cudaFree(_result_dev);
  cudaFree(_adjacency_dev);
  if(_num_degrees_flag)
    cudaFree(_num_degrees_dev);

  if(_rank == 0 && _is_profile)
    ODP_Profile("MPI+CUDA", _kind, _symmetries, _mem_usage,
		_elapsed_time, _times, _procs);
}

extern "C" void ODP_Set_aspl_mpi_cuda(const int* __restrict__ adjacency,
				      int *diameter, long *sum, double *ASPL)
{
  double t = ODP_Get_time();
  
  cudaMemcpy(_adjacency_dev, adjacency, sizeof(int)*(_nodes/_symmetries)*_degree, cudaMemcpyHostToDevice);
  
  if(_kind == ASPL_NORMAL)
    aspl_mpi_cuda_mat       (adjacency, diameter, sum, ASPL);
  else
    aspl_mpi_cuda_mat_saving(adjacency, diameter, sum, ASPL);

  _elapsed_time += ODP_Get_time() - t;
    
  if(*diameter > _nodes)
    ERROR("This graph is not connected graph.\n");

  _times++;
}
