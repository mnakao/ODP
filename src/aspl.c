#include "common.h"
static uint64_t *_A, *_B;
static int _nodes, _degree, _symmetries, _kind, _width, _height;
static int* _num_degrees = NULL, *_itable = NULL;
static int* _frontier = NULL, *_distance = NULL, *_next = NULL;
static char* _bitmap = NULL;
static unsigned int _elements, _times;
static double _mem_usage, _elapsed_time;
static bool _enable_avx2 = false, _is_profile = false, _enable_grid_s = false;

extern void ODP_Create_itable(const int width, const int height, const int symmetries, int *_itable);
extern bool ODP_Check_profile();
extern double ODP_Get_time();
extern int ODP_LOCAL_INDEX_GRID(const int x, const int width, const int height, const int symmetries);
extern int ODP_ROTATE(const int v, const int width, const int height, const int symmetries, const int degree);
extern void ODP_Profile(const char* name, const int kind, const int symmetries, const double mem_usage,
			const double elapsed_time, const unsigned int times, const int procs);
extern int ODP_Get_kind(const int nodes, const int degree, const int* num_degrees, const int symmetries,
			const int procs, const bool is_cpu);
extern double ODP_Get_mem_usage(const int kind, const int nodes, const int degree, const int symmetries,
				const int *num_degrees, const int procs, const bool is_cpu);
extern void ODP_Matmul(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int height, const int degree,
		       const int *restrict num_degrees, const int *restrict adjacency, const bool enable_avx2, const int elements, const int symmetries, const int itable[nodes]);
extern void ODP_Matmul_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int height, const int degree,
			     const int *restrict num_degrees, const int *restrict adjacency, const bool enable_avx2, const int symmetries, const int itable[nodes]);
extern void ODP_Malloc(uint64_t **a, const size_t s, const bool enable_avx2);
extern void ODP_Free(uint64_t *a, const bool enable_avx2);
extern int ODP_top_down_step(const int level, const int num_frontier, const int* restrict adjacency,
			     const int nodes, const int degree, const int* restrict num_degrees, const bool enable_grid_s,
			     const int width, const int height, const int symmetries,
			     int* restrict frontier, int* restrict next, int* restrict distance, char* restrict bitmap);
  
static void aspl_mat(const int* restrict adjacency,
		     int *diameter, long *sum, double *ASPL)
{
#pragma omp parallel for
  for(int i=0;i<_nodes*_elements;i++)
    _A[i] = _B[i] = 0;

#pragma omp parallel for
  for(int i=0;i<_nodes/_symmetries;i++){
    unsigned int offset = i*_elements+i/UINT64_BITS;
    _A[offset] = _B[offset] = (0x1ULL << (i%UINT64_BITS));
  }

  *sum = (long)_nodes * (_nodes - 1);
  *diameter = 1;
  for(int kk=0;kk<_nodes;kk++){
    ODP_Matmul(_A, _B, _nodes, _height, _degree, _num_degrees, adjacency, 
	       _enable_avx2, _elements, _symmetries, _itable);

    uint64_t num = 0;
#pragma omp parallel for reduction(+:num)
    for(int i=0;i<_elements*_nodes;i++)
      num += POPCNT(_B[i]);

    num *= _symmetries;
    if(num == (uint64_t)_nodes*_nodes) break;

    // swap A <-> B
    uint64_t* tmp = _A;
    _A = _B;
    _B = tmp;

    *sum += (long)_nodes * _nodes - num;
    (*diameter) += 1;
  }

  *ASPL = *sum / (((double)_nodes-1)*_nodes);
  *sum /= 2.0;
}

static void aspl_mat_saving(const int* restrict adjacency,
			    int *diameter, long *sum, double *ASPL)
{
  int parsize = (_elements+(CPU_CHUNK-1))/CPU_CHUNK;

  *sum = (long)_nodes * (_nodes - 1);
  *diameter = 1;
  for(int t=0;t<parsize;t++){
    unsigned int kk, l;
#pragma omp parallel for
    for(int i=0;i<_nodes*CPU_CHUNK;i++)
      _A[i] = _B[i] = 0;
    
    for(l=0; l<UINT64_BITS*CPU_CHUNK && UINT64_BITS*t*CPU_CHUNK+l<_nodes/_symmetries; l++){
      unsigned int offset = (UINT64_BITS*t*CPU_CHUNK+l)*CPU_CHUNK+l/UINT64_BITS;
      _A[offset] = _B[offset] = (0x1ULL<<(l%UINT64_BITS));
    }

    for(kk=0;kk<_nodes;kk++){
      ODP_Matmul_CHUNK(_A, _B, _nodes, _height, _degree, _num_degrees, adjacency, _enable_avx2, _symmetries, _itable);

      uint64_t num = 0;
#pragma omp parallel for reduction(+:num)
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

  *ASPL = *sum / (((double)_nodes-1)*_nodes);
  *sum /= 2.0;
}

static void init_aspl_s(const int nodes, const int degree, const int* num_degrees, const int symmetries)
{
  if(nodes % symmetries != 0)
    ERROR("nodes(%d) must be divisible by symmetries(%d)\n", nodes, symmetries);
  else if(CPU_CHUNK % 4 != 0)
    ERROR("CPU_CHUNK(%d) in parameter.h must be multiple of 4\n", CPU_CHUNK);

  _kind      = ODP_Get_kind(nodes, degree, num_degrees, symmetries, 1, true);
  _mem_usage = ODP_Get_mem_usage(_kind, nodes, degree, symmetries, num_degrees, 1, true);
  _elements  = (nodes/symmetries+(UINT64_BITS-1))/UINT64_BITS;
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
  else{ // _kind == ASPL_BFS
    _bitmap   = malloc(sizeof(char) * nodes);
    _frontier = malloc(sizeof(int)  * nodes);
    _distance = malloc(sizeof(int)  * nodes);
    _next     = malloc(sizeof(int)  * nodes);
  }

  _nodes = nodes;
  _degree = degree;
  _symmetries = symmetries;
  _is_profile = ODP_Check_profile();
  _elapsed_time = 0;
  _times = 0;

  if(num_degrees){
    _num_degrees = malloc(sizeof(int) * nodes);
    memcpy(_num_degrees, num_degrees, sizeof(int) * nodes);
  }
}

static void aspl_bfs(const int* restrict adjacency, int* diameter, long *sum, double* ASPL)
{
  int based_nodes = _nodes/_symmetries;
  bool reached = true;
  *diameter = 0;
  *sum      = 0;

  for(int s=0;s<based_nodes;s++){
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
      _frontier[0]  = s;
      _distance[s] = level;
      _bitmap[s]   = VISITED;
    }

    while(1){
      num_frontier = ODP_top_down_step(level++, num_frontier, adjacency, _nodes, _degree, _num_degrees,
				       _enable_grid_s, _width, _height, _symmetries, _frontier, _next, _distance, _bitmap);
      if(num_frontier == 0) break;

      int *tmp = _frontier;
      _frontier = _next;
      _next     = tmp;
    }

    *diameter = MAX(*diameter, level-1);

    if(s == 0){
      for(int i=1;i<_nodes;i++)
        if(_bitmap[i] == NOT_VISITED)
          reached = false;
      
      if(!reached){
        *diameter = INT_MAX;
        return;
      }
    }

    for(int i=0;i<_nodes;i++)
      *sum += (_distance[i] + 1) * _symmetries;
  }
  
  *sum = (*sum - _nodes)/2;
  *ASPL = *sum / (((double)_nodes-1)*_nodes) * 2;
}

void ODP_Init_aspl_general(const int nodes, const int degree, const int* num_degrees)
{
  init_aspl_s(nodes, degree, num_degrees, 1);
}

void ODP_Init_aspl_general_s(const int nodes, const int degree, const int* num_degrees, const int symmetries)
{

  if(num_degrees){
    int *tmp_num_degrees = malloc(sizeof(int) * nodes);
    int based_nodes = nodes/symmetries;
    for(int i=0;i<symmetries;i++)
      for(int j=0;j<based_nodes;j++)
        tmp_num_degrees[i*based_nodes+j] = num_degrees[j];

    init_aspl_s(nodes, degree, tmp_num_degrees, symmetries);
    free(tmp_num_degrees);
  }
  else{
    init_aspl_s(nodes, degree, NULL, symmetries);
  }
}

void ODP_Init_aspl_grid(const int width, const int height, const int degree, const int* num_degrees)
{
  int nodes = width * height;
  _width  = width;
  _height = height;
  init_aspl_s(nodes, degree, num_degrees, 1);
}

void ODP_Init_aspl_grid_s(const int width, const int height, const int degree, const int* num_degrees, const int symmetries)

{
  int nodes = width * height;
  _width  = width;
  _height = height;
  if(symmetries == 2 || symmetries == 4)
    _enable_grid_s = true;
  
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
    init_aspl_s(nodes, degree, tmp_num_degrees, symmetries);
    free(tmp_num_degrees);
  }
  else{
    init_aspl_s(nodes, degree, NULL, symmetries);
  }

  if(symmetries > 1){
    _itable = malloc(sizeof(int) * nodes);
    ODP_Create_itable(width, height, symmetries, _itable);
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
  }
  
  if(_is_profile){
#ifdef _OPENMP
    ODP_Profile("THREADS", _kind, _symmetries, _mem_usage,
		 _elapsed_time, _times, 1);
#else
    ODP_Profile("SERIAL",  _kind, _symmetries, _mem_usage,
		 _elapsed_time, _times, 1);
#endif
  }
}

void ODP_Set_aspl(const int* restrict adjacency, int *diameter, long *sum, double *ASPL)
{
  double t = ODP_Get_time();

  if(_kind == ASPL_MATRIX)
    aspl_mat       (adjacency, diameter, sum, ASPL);
  else if(_kind == ASPL_MATRIX_SAVING)
    aspl_mat_saving(adjacency, diameter, sum, ASPL);
  else // _kind == ASPL_MATRIX_BFS
    aspl_bfs(adjacency, diameter, sum, ASPL);

  _elapsed_time += ODP_Get_time() - t;
  
  if(*diameter > _nodes){
    *diameter = INT_MAX;
    *sum = LONG_MAX;
    *ASPL = DBL_MAX;
  }
  
  _times++;
}
