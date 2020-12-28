#include "common.h"
#include <mpi.h>
static uint64_t *_A, *_B;
static int _nodes, _degree, _groups, _kind, _rank, _procs;
static const int* _num_degrees;
static double _mem_usage;
static MPI_Comm _comm;

extern int  apsp_get_kind(const int nodes, const int degree, const int* num_degrees, const int groups,
                          const int procs);
extern void apsp_start_profile();
extern void apsp_end_profile(const char* name, const int kind, const int groups, const double mem_usage, const int procs);
extern double apsp_get_mem_usage(const int kind, const int nodes, const int degree, const int groups,
				 const int *num_degree, const int procs);
  
static void apsp_mpi_mat(const int* restrict adjacency,
			 int *diameter, long *sum, double *ASPL)
{
  unsigned int elements = (_nodes/_groups+(UINT64_BITS-1))/UINT64_BITS;
  unsigned int chunk = (elements+(_procs-1))/_procs;
  int parsize = (elements+(chunk-1))/chunk;

  *sum = 0.0;
  *diameter = 1;
  for(int t=_rank;t<parsize;t+=_procs){
    uint64_t kk, l;
#pragma omp parallel for
    for(int i=0;i<_nodes*chunk;i++)
      _A[i] = _B[i] = 0;

    for(l=0; l<UINT64_BITS*chunk && UINT64_BITS*t*chunk+l<_nodes/_groups; l++){
      unsigned int offset = (UINT64_BITS*t*chunk+l)*chunk+l/UINT64_BITS;
      _A[offset] = _B[offset] = (0x1ULL<<(l%UINT64_BITS));
    }

    for(kk=0;kk<_nodes;kk++){
      if(!_num_degrees){
#pragma omp parallel for
	for(int i=0;i<_nodes;i++)
	  for(int j=0;j<_degree;j++){
	    int n = *(adjacency + i * _degree + j);  // int n = adjacency[i][j];
	    for(int k=0;k<chunk;k++)
	      _B[i*chunk+k] |= _A[n*chunk+k];
	  }
      }
      else{
#pragma omp parallel for
	for(int i=0;i<_nodes;i++)
	  for(int j=0;j<_num_degrees[i];j++){
	    int n = *(adjacency + i * _degree + j);  // int n = adjacency[i][j];
	    for(int k=0;k<chunk;k++)
	      _B[i*chunk+k] |= _A[n*chunk+k];
	  }
      }

      uint64_t num = 0;
#pragma omp parallel for reduction(+:num)
      for(int i=0;i<chunk*_nodes;i++)
        num += POPCNT(_B[i]);

      if(num == (uint64_t)_nodes*l) break;

      // swap A <-> B
      uint64_t* tmp = _A;
      _A = _B;
      _B = tmp;

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

static void apsp_mpi_mat_saving(const int* restrict adjacency,
				int *diameter, long *sum, double *ASPL)
{
  unsigned int elements = (_nodes/_groups+(UINT64_BITS-1))/UINT64_BITS;
  int parsize = (elements+(CHUNK-1))/CHUNK;
 
  *sum = 0.0;
  *diameter = 1;
  for(int t=_rank;t<parsize;t+=_procs){
    unsigned int kk, l;
#pragma omp parallel for
    for(int i=0;i<_nodes*CHUNK;i++)
      _A[i] = _B[i] = 0;
    
    for(l=0; l<UINT64_BITS*CHUNK && UINT64_BITS*t*CHUNK+l<_nodes/_groups; l++){
      unsigned int offset = (UINT64_BITS*t*CHUNK+l)*CHUNK+l/UINT64_BITS;
      _A[offset] = _B[offset] = (0x1ULL<<(l%UINT64_BITS));
    }

    for(kk=0;kk<_nodes;kk++){
      if(!_num_degrees){
#pragma omp parallel for
	for(int i=0;i<_nodes;i++)
	  for(int j=0;j<_degree;j++){
	    int n = *(adjacency + i * _degree + j);  // int n = adjacency[i][j];
	    for(int k=0;k<CHUNK;k++)
	      _B[i*CHUNK+k] |= _A[n*CHUNK+k];
	  }
      }
      else{
#pragma omp parallel for
        for(int i=0;i<_nodes;i++)
      	  for(int j=0;j<_num_degrees[i];j++){
            int n = *(adjacency + i * _degree + j);  // int n = adjacency[i][j];
            for(int k=0;k<CHUNK;k++)
              _B[i*CHUNK+k] |= _A[n*CHUNK+k];
          }
      }

      uint64_t num = 0;
#pragma omp parallel for reduction(+:num)
      for(int i=0;i<CHUNK*_nodes;i++)
        num += POPCNT(_B[i]);

      if(num == (uint64_t)_nodes*l) break;
      
      // swap A <-> B
      uint64_t* tmp = _A;
      _A = _B;
      _B = tmp;

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

void apsp_mpi_init_s(const int nodes, const int degree,
		     const int* restrict num_degrees, const MPI_Comm comm, const int groups)
{
  if(nodes % groups != 0)
    ERROR("nodes(%d) must be divisible by group(%d)\n", nodes, groups);

  MPI_Comm_rank(comm, &_rank);
  MPI_Comm_size(comm, &_procs);

  _kind = apsp_get_kind(nodes, degree, num_degrees, groups, _procs);
  _mem_usage = apsp_get_mem_usage(_kind, nodes, degree, groups, num_degrees, _procs);
  
  if(_kind == APSP_NORMAL){
    unsigned int elements = (nodes/groups+(UINT64_BITS-1))/UINT64_BITS;
    unsigned int chunk = (elements+(_procs-1))/_procs;
    size_t s = nodes * chunk * sizeof(uint64_t);
    posix_memalign((void **)&_A, ALIGN_VALUE, s); // uint64_t A[nodes][chunk];
    posix_memalign((void **)&_B, ALIGN_VALUE, s); // uint64_t B[nodes][chunk];
  }
  else{
    size_t s = nodes * CHUNK * sizeof(uint64_t);
    posix_memalign((void **)&_A, ALIGN_VALUE, s); // uint64_t A[nodes][CHUNK];
    posix_memalign((void **)&_B, ALIGN_VALUE, s); // uint64_t B[nodes][CHUNK];
  }
  _nodes = nodes;
  _degree = degree;
  _num_degrees = num_degrees;
  _groups = groups;
  _comm = comm;
}

void apsp_mpi_init(const int nodes, const int degree,
		   const int* restrict num_degrees, MPI_Comm comm)
{
  apsp_mpi_init_s(nodes, degree, num_degrees, comm, 1);
}

void apsp_mpi_finalize()
{
  free(_A);
  free(_B);
}

void apsp_mpi_run(const int* restrict adjacency, int *diameter, long *sum, double *ASPL)
{
  MPI_Barrier(_comm);
  if(_rank == 0)
    apsp_start_profile();
  
  if(_kind == APSP_NORMAL)
    apsp_mpi_mat       (adjacency, diameter, sum, ASPL);
  else
    apsp_mpi_mat_saving(adjacency, diameter, sum, ASPL);

  if(*diameter > _nodes)
    ERROR("This graph is not connected graph.\n");

  if(_rank == 0){
#ifdef _OPENMP
    apsp_end_profile("MPI+Threads", _kind, _groups, _mem_usage, _procs);
#else
    apsp_end_profile("MPI", _kind, _groups, _mem_usage, _procs);
#endif
  }
}
