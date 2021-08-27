#include "common.h"
#include <mpi.h>
static uint64_t *_A, *_B;
static int _nodes, _degree, _symmetries, _kind, _rank, _procs, _height = -1;
static int* _num_degrees = NULL, *_itable = NULL;
static int* _frontier = NULL, *_distance = NULL, *_next = NULL;
static char* _bitmap = NULL;
static unsigned int _elements, _total_elements, _times;
static double _mem_usage, _elapsed_time;
static bool _enable_avx2 = false, _is_profile = false, _enable_grid_s = false;
static MPI_Comm _comm;

extern bool ODP_Check_profile();
extern double ODP_Get_time();
extern void ODP_Create_itable(const int width, const int height, const int symmetries, int *_itable);
extern int ODP_LOCAL_INDEX_GRID(const int x, const int width, const int height, const int symmetries);
extern int ODP_ROTATE(const int v, const int width, const int height, const int symmetries, const int degree);
extern void ODP_Profile(const char* name, const int kind, const int symmetries, const double mem_usage,
			const double elapsed_time, const unsigned int times, const int procs);
extern int ODP_Get_kind(const int nodes, const int degree, const int* num_degrees, const int symmetries,
			const int procs, const bool is_cpu, const bool enable_grid_s);
extern double ODP_Get_mem_usage(const int kind, const int nodes, const int degree, const int symmetries,
				const int *num_degrees, const int procs, const bool is_cpu, const bool enable_grid_s);
extern void ODP_Matmul(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int height,
		       const int degree, const int *restrict num_degrees, const int *restrict adjacency,
		       const int *itable, const int elements, const int symmetries, const bool enable_avx2);
extern void ODP_Matmul_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int height,
			     const int degree, const int *num_degrees, const int *restrict adjacency,
			     const int *itable, const int symmetries, const bool enable_avx2);
extern void ODP_Malloc(uint64_t **a, const size_t s, const bool enable_avx2);
extern void ODP_Free(uint64_t *a, const bool enable_avx2);
extern int ODP_top_down_step(const int level, const int num_frontier, const int* restrict adjacency,
                             const int nodes, const int degree, const int* restrict num_degrees, const bool enable_grid_s,
                             const int height, const int symmetries,
                             int* restrict frontier, int* restrict next, int* restrict distance, char* restrict bitmap);

static void aspl_mpi_mat(const int* restrict adjacency,
			 int *diameter, long *sum, double *ASPL)
{
  int parsize = (_total_elements+(_elements-1))/_elements;

  *sum = 0.0;
  *diameter = 1;
  for(int t=_rank;t<parsize;t+=_procs){
    uint64_t kk, l;
#pragma omp parallel for
    for(int i=0;i<_nodes*_elements;i++)
      _A[i] = _B[i] = 0;

    for(l=0; l<UINT64_BITS*_elements && UINT64_BITS*t*_elements+l<_nodes/_symmetries; l++){
      unsigned int offset = (UINT64_BITS*t*_elements+l)*_elements+l/UINT64_BITS;
      _A[offset] = _B[offset] = (0x1ULL<<(l%UINT64_BITS));
    }

    for(kk=0;kk<_nodes;kk++){
      ODP_Matmul(_A, _B, _nodes, _height, _degree, _num_degrees, adjacency,
		 _itable, _elements, _symmetries, _enable_avx2);

      uint64_t num = 0;
#ifndef __FUGAKU
#pragma omp parallel for reduction(+:num)
#endif
      for(int i=0;i<_elements*_nodes;i++)
        num += POPCNT(_B[i]);

      if(num == (uint64_t)_nodes*l) break;

      // swap A <-> B
      uint64_t* tmp = _A;
      _A = _B;
      _B = tmp;

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

static void aspl_mpi_mat_saving(const int* restrict adjacency,
				int *diameter, long *sum, double *ASPL)
{
  int parsize = (_total_elements+(CPU_CHUNK-1))/CPU_CHUNK;
 
  *sum = 0.0;
  *diameter = 1;
  for(int t=_rank;t<parsize;t+=_procs){
    unsigned int kk, l;
#pragma omp parallel for
    for(int i=0;i<_nodes*CPU_CHUNK;i++)
      _A[i] = _B[i] = 0;

    for(l=0; l<UINT64_BITS*CPU_CHUNK && UINT64_BITS*t*CPU_CHUNK+l<_nodes/_symmetries; l++){
      unsigned int offset = (UINT64_BITS*t*CPU_CHUNK+l)*CPU_CHUNK+l/UINT64_BITS;
      _A[offset] = _B[offset] = (0x1ULL<<(l%UINT64_BITS));
    }

    for(kk=0;kk<_nodes;kk++){
      ODP_Matmul_CHUNK(_A, _B, _nodes, _height, _degree, _num_degrees, adjacency, _itable,
		       _symmetries, _enable_avx2);

      uint64_t num = 0;
#ifndef __FUGAKU
#pragma omp parallel for reduction(+:num)
#endif
      for(int i=0;i<CPU_CHUNK*_nodes;i++)
        num += POPCNT(_B[i]);

      if(num == (uint64_t)_nodes*l) break;
      
      // swap A <-> B
      uint64_t* tmp = _A;
      _A = _B;
      _B = tmp;

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

static void init_aspl_mpi_s(const int nodes, const int degree,
			    const int* restrict num_degrees, const MPI_Comm comm, const int symmetries)
{
  if(nodes % symmetries != 0)
    ERROR("nodes(%d) must be divisible by symmetries(%d)\n", nodes, symmetries);
  else if(CPU_CHUNK % 4 != 0)
    ERROR("CPU_CHUNK(%d) in parameter.h must be multiple of 4\n", CPU_CHUNK);

  MPI_Comm_rank(comm, &_rank);
  MPI_Comm_size(comm, &_procs);

  _kind = ODP_Get_kind(nodes, degree, num_degrees, symmetries, _procs, true, _enable_grid_s);
  _mem_usage = ODP_Get_mem_usage(_kind, nodes, degree, symmetries, num_degrees, _procs, true, _enable_grid_s);
  _total_elements = (nodes/symmetries+(UINT64_BITS-1))/UINT64_BITS;
  _elements = (_total_elements+(_procs-1))/_procs;
#ifdef __AVX2__
  if(_elements >= 4){ // For performance
    _enable_avx2 = true;
    _elements = ((_elements+3)/4)*4;  // _elements must be multiple of 4
  }
#endif

  if(_kind == ASPL_MATRIX || _kind == ASPL_MATRIX_SAVING){
    size_t s = (_kind == ASPL_MATRIX)? _elements : CPU_CHUNK;
    ODP_Malloc(&_A, nodes*s*sizeof(uint64_t), _enable_avx2); // uint64_t A[nodes][s];
    ODP_Malloc(&_B, nodes*s*sizeof(uint64_t), _enable_avx2); // uint64_t B[nodes][s];
  }
  else{
    _bitmap   = malloc(sizeof(char) * nodes);
    _frontier = malloc(sizeof(int)  * nodes);
    _distance = malloc(sizeof(int)  * nodes);
    _next     = malloc(sizeof(int)  * nodes);
#ifdef _OPENMP
    ODP_declare_local_frontier(nodes);
#endif
  }
  
  _nodes = nodes;
  _degree = degree;
  _symmetries = symmetries;
  _comm = comm;
  _is_profile = ODP_Check_profile();
  _elapsed_time = 0;
  _times = 0;

  if(num_degrees){
    _num_degrees = malloc(sizeof(int) * nodes);
    memcpy(_num_degrees, num_degrees, sizeof(int) * nodes);
  }
}

static void aspl_mpi_bfs(const int* restrict adjacency, int* diameter, long *sum, double* ASPL)
{
  int based_nodes = _nodes/_symmetries;
  bool reached = true;
  *diameter = 0;
  *sum      = 0;

  for(int s=_rank;s<based_nodes;s+=_procs){
    int num_frontier = 1, level = 0;
    for(int i=0;i<_nodes;i++)
      _bitmap[i] = NOT_VISITED;

    if(_enable_grid_s && _symmetries == 4){
      int based_height = _height/2;
      int ss = (s/based_height) * _height + (s%based_height);
      _frontier[0]  = ss;
      _distance[ss] = level;
      _bitmap[ss]   = VISITED;
    }
    else{
      _frontier[0] = s;
      _distance[s] = level;
      _bitmap[s]   = VISITED;
    }

    while(1){
      num_frontier = ODP_top_down_step(level++, num_frontier, adjacency, _nodes, _degree, _num_degrees,
                                       _enable_grid_s, _height, _symmetries, _frontier, _next, _distance, _bitmap);
      if(num_frontier == 0) break;

      int *tmp = _frontier;
      _frontier = _next;
      _next     = tmp;
    }

    *diameter = MAX(*diameter, level-1);

    for(int i=0;i<_nodes;i++)
      *sum += (_distance[i] + 1) * _symmetries;
  }

  MPI_Bcast(&reached, 1, MPI_C_BOOL, 0, _comm);
  if(!reached){
    *diameter = INT_MAX;
    return;
  }
  MPI_Allreduce(MPI_IN_PLACE, diameter, 1, MPI_INT,  MPI_MAX, _comm);
  MPI_Allreduce(MPI_IN_PLACE, sum,      1, MPI_LONG, MPI_SUM, _comm);
  *sum = (*sum - _nodes)/2;
  *ASPL = *sum / (((double)_nodes-1)*_nodes) * 2;
}

void ODP_Init_aspl_mpi_general(const int nodes, const int degree,
			       const int* restrict num_degrees, const MPI_Comm comm)
{
  init_aspl_mpi_s(nodes, degree, num_degrees, comm, 1);
}

void ODP_Init_aspl_mpi_general_s(const int nodes, const int degree,
				 const int* restrict num_degrees, const MPI_Comm comm, const int symmetries)
{
  if(num_degrees){
    int *tmp_num_degrees = malloc(sizeof(int) * nodes);
    int based_nodes = nodes/symmetries;
    for(int i=0;i<symmetries;i++)
      for(int j=0;j<based_nodes;j++)
        tmp_num_degrees[i*based_nodes+j] = num_degrees[j];
    
    init_aspl_mpi_s(nodes, degree, tmp_num_degrees, comm, symmetries);
    free(tmp_num_degrees);
  }
  else{
    init_aspl_mpi_s(nodes, degree, NULL, comm, symmetries);
  }
}

void ODP_Init_aspl_mpi_grid(const int width, const int height, const int degree,
			    const int* restrict num_degrees, const MPI_Comm comm)
{
  int nodes = width * height;
  _height = height;
  init_aspl_mpi_s(nodes, degree, num_degrees, comm, 1);
}

void ODP_Init_aspl_mpi_grid_s(const int width, const int height, const int degree, const int* num_degrees,
			      const MPI_Comm comm, const int symmetries)
{
  int nodes = width * height;
  _height = height;
  if(symmetries == 2 || symmetries == 4){
    _enable_grid_s = true;
    _itable = malloc(sizeof(int) * nodes);
    ODP_Create_itable(width, height, symmetries, _itable);
  }
  
  if(num_degrees){
    int *tmp_num_degrees = malloc(sizeof(int) * nodes);
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
    init_aspl_mpi_s(nodes, degree, tmp_num_degrees, comm, symmetries);
    free(tmp_num_degrees);
  }
  else{
    init_aspl_mpi_s(nodes, degree, NULL, comm, symmetries);
  }
}

void ODP_Finalize_aspl()
{
  if(_kind == ASPL_MATRIX || _kind == ASPL_MATRIX_SAVING){
    ODP_Free(_A, _enable_avx2);
    ODP_Free(_B, _enable_avx2);
    if(_num_degrees) free(_num_degrees);
    if(_itable)      free(_itable);
  }
  else{ // _kind == ASPL_BFS
    free(_bitmap);
    free(_frontier);
    free(_distance);
    free(_next);
#ifdef _OPENMP
    ODP_free_local_frontier();
#endif
  }
  
  if(_rank == 0 && _is_profile){
#ifdef _OPENMP
    ODP_Profile("MPI+THREADS", _kind, _symmetries, _mem_usage,
		_elapsed_time, _times, _procs);
#else
    ODP_Profile("MPI", _kind, _symmetries, _mem_usage,
		_elapsed_time, _times, _procs);
#endif
  }
}

void ODP_Set_aspl(const int* restrict adjacency, int *diameter, long *sum, double *ASPL)
{
  double t = ODP_Get_time();
  
  if(_kind == ASPL_MATRIX)
    aspl_mpi_mat       (adjacency, diameter, sum, ASPL);
  else if(_kind == ASPL_MATRIX_SAVING)
    aspl_mpi_mat_saving(adjacency, diameter, sum, ASPL);
  else // _kind == ASPL_MATRIX_BFS
    aspl_mpi_bfs(adjacency, diameter, sum, ASPL);
  
  _elapsed_time += ODP_Get_time() - t;
    
  if(*diameter > _nodes){
    *diameter = INT_MAX;
    *sum = LONG_MAX;
    *ASPL = DBL_MAX;
  }

  _times++;
}
