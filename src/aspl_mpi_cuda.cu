#include "common.h"
#include <mpi.h>
static uint64_t *_A_dev, *_B_dev;
static uint64_t *_result, *_result_dev;
static int *_adjacency_dev, *_num_degrees_dev = NULL;
static bool _is_profile = false, _enable_grid_s = false;;
static int _nodes, _degree, _symmetries, _rank, _procs, _kind, _height = -1;
static double _mem_usage, _elapsed_time;
static unsigned int _times;
static MPI_Comm _comm;

extern "C" bool ODP_Check_profile();
extern "C" double ODP_Get_time();
extern "C" int ODP_LOCAL_INDEX_GRID(const int x, const int width, const int height, const int symmetries);
extern "C" int ODP_ROTATE(const int v, const int width, const int height, const int symmetries, const int degree);
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
				       const int* __restrict__ num_degrees, const int nodes, const int height, const int degree,
				       const unsigned int elements, const int symmetries, const int enable_grid_s);
extern __global__ void ODP_Matmul_CHUNK_cuda(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
                                             const int* __restrict__ num_degrees, const int nodes, const int _height, const int degree,
                                             const int symmetries, const bool enable_grid_s);

static __global__ void init_buffers(uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int nodes,
				    const int symmetries, const int t, const int chunk, const int height, const bool enable_grid_s)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(enable_grid_s && symmetries == 4){
    int based_height = height/2;
    while (tid<UINT64_BITS*chunk && UINT64_BITS*t*chunk+tid<nodes/symmetries) {
      int i = (tid/based_height) * height + (tid%based_height);
      unsigned int offset = (UINT64_BITS*t*chunk+i)*chunk+i/UINT64_BITS;
      A[offset] = B[offset] = (0x1ULL<<(i%UINT64_BITS));
      tid += blockDim.x * gridDim.x;
    }
  }
  else{
    while (tid<UINT64_BITS*chunk && UINT64_BITS*t*chunk+tid<nodes/symmetries) {
      unsigned int offset = (UINT64_BITS*t*chunk+tid)*chunk+tid/UINT64_BITS;
      A[offset] = B[offset] = (0x1ULL<<(tid%UINT64_BITS));
      tid += blockDim.x * gridDim.x;
    }
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
    init_buffers      <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes, _symmetries, t, chunk, _height, _enable_grid_s);

    for(kk=0;kk<_nodes;kk++){
      ODP_Matmul_cuda <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _adjacency_dev, _num_degrees_dev,
					       _nodes, _height, _degree, chunk, _symmetries, _enable_grid_s);
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
    init_buffers      <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _nodes, _symmetries, t, GPU_CHUNK, _height, _enable_grid_s);

    for(kk=0;kk<_nodes;kk++){
      ODP_Matmul_CHUNK_cuda <<< BLOCKS, THREADS >>> (_A_dev, _B_dev, _adjacency_dev, _num_degrees_dev,
						     _nodes, _height, _degree, _symmetries, _enable_grid_s);
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

static void init_aspl_mpi_cuda_s(const int nodes, const int degree, const int* __restrict__ num_degrees,
				 const MPI_Comm comm, const int symmetries)
{
  cuInit(0);
  
  if(nodes % symmetries != 0)
    ERROR("nodes(%d) must be divisible by symmetries(%d)\n", nodes, symmetries);

  MPI_Comm_rank(comm, &_rank);
  MPI_Comm_size(comm, &_procs);
  _kind = ODP_Get_kind(nodes, degree, num_degrees, symmetries, _procs, false);
  _mem_usage = ODP_Get_mem_usage(_kind, nodes, degree, symmetries, num_degrees, _procs, false);
  size_t elements = (nodes/symmetries+(UINT64_BITS-1))/UINT64_BITS;
  size_t s = (_kind == ASPL_MATRIX)? (elements+(_procs-1))/_procs : GPU_CHUNK;
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
  }
  _is_profile = ODP_Check_profile();
  _elapsed_time = 0;
  _times = 0;
}

extern "C" void ODP_Init_aspl_mpi_cuda_general(const int nodes, const int degree,
					       const int* __restrict__ num_degrees, const MPI_Comm comm)
{
  init_aspl_mpi_cuda_s(nodes, degree, num_degrees, comm, 1);
}

extern "C" void ODP_Init_aspl_mpi_cuda_general_s(const int nodes, const int degree, const int* __restrict__ num_degrees,
						 const MPI_Comm comm, const int symmetries)
{
  if(num_degrees){
    int *tmp_num_degrees = (int *)malloc(sizeof(int) * nodes);
    int based_nodes = nodes/symmetries;
    for(int i=0;i<symmetries;i++)
      for(int j=0;j<based_nodes;j++)
        tmp_num_degrees[i*based_nodes+j] = num_degrees[j];

    init_aspl_mpi_cuda_s(nodes, degree, tmp_num_degrees, comm, symmetries);
    free(tmp_num_degrees);
  }
  else{
    init_aspl_mpi_cuda_s(nodes, degree, NULL, comm, symmetries);
  }
}

extern "C" void ODP_Init_aspl_mpi_cuda_grid(const int width, const int height, const int degree,
					    const int* __restrict__ num_degrees, const MPI_Comm comm)
{
  int nodes = width * height;
  _height = height;
  init_aspl_mpi_cuda_s(nodes, degree, num_degrees, comm, 1);
}

extern "C" void ODP_Init_aspl_mpi_cuda_grid_s(const int width, const int height, const int degree,
					      const int* __restrict__ num_degrees, const MPI_Comm comm, const int symmetries)
{
  int nodes = width * height;
  _height = height;
  if(symmetries == 2 || symmetries == 4)
    _enable_grid_s = true;
    
  if(num_degrees){
    int *tmp_num_degrees = (int *)malloc(sizeof(int) * nodes);
    int based_nodes = nodes/symmetries;
    if(symmetries == 2){
      for(int i=0;i<based_nodes;i++){
        tmp_num_degrees[i] = num_degrees[i];
        tmp_num_degrees[ODP_ROTATE(i, width, height, symmetries, 180)] = num_degrees[i];
      }
    }
    else if(symmetries == 4){
      for(int i=0;i<based_nodes;i++){
        int v = ODP_LOCAL_INDEX_GRID(i,width,height,symmetries);
        tmp_num_degrees[v] = num_degrees[i];
        tmp_num_degrees[ODP_ROTATE(v, width, height, symmetries,  90)] = num_degrees[i];
        tmp_num_degrees[ODP_ROTATE(v, width, height, symmetries, 180)] = num_degrees[i];
        tmp_num_degrees[ODP_ROTATE(v, width, height, symmetries, 270)] = num_degrees[i];
      }
    }
    init_aspl_mpi_cuda_s(nodes, degree, tmp_num_degrees, comm, symmetries);
    free(tmp_num_degrees);
  }
  else{
    init_aspl_mpi_cuda_s(nodes, degree, NULL, comm, symmetries);
  }
}

extern "C" void ODP_Finalize_aspl()
{
  cudaFree(_A_dev);
  cudaFree(_B_dev);
  cudaFreeHost(_result);
  cudaFree(_result_dev);
  cudaFree(_adjacency_dev);
  if(_num_degrees_dev) cudaFree(_num_degrees_dev);

  if(_rank == 0 && _is_profile)
    ODP_Profile("MPI+CUDA", _kind, _symmetries, _mem_usage,
		_elapsed_time, _times, _procs);
}

extern "C" void ODP_Set_aspl(const int* __restrict__ adjacency, int *diameter, long *sum, double *ASPL)
{
  double t = ODP_Get_time();
  
  cudaMemcpy(_adjacency_dev, adjacency, sizeof(int)*(_nodes/_symmetries)*_degree, cudaMemcpyHostToDevice);
  
  if(_kind == ASPL_MATRIX)
    aspl_mpi_cuda_mat       (adjacency, diameter, sum, ASPL);
  else
    aspl_mpi_cuda_mat_saving(adjacency, diameter, sum, ASPL);

  _elapsed_time += ODP_Get_time() - t;
    
  if(*diameter > _nodes){
    *diameter = INT_MAX;
    *sum = LONG_MAX;
    *ASPL = DBL_MAX;
  }

  _times++;
}
