#include "common.h"
static uint64_t *_A, *_B;
static int _nodes, _degree, _symmetries, _kind;
static const int* _num_degrees;
static unsigned int _elements, _times;
static double _mem_usage, _elapsed_time;
static bool _enable_avx2 = false, _is_profile;
static time_t _start_t;

extern bool apsp_check_profile();
extern double apsp_get_time();
extern void apsp_profile(const char* name, const int kind, const int symmetries, const double mem_usage, const time_t start_t,
                         const time_t end_t, const double elapsed_time, const unsigned int times, const int procs);
extern int apsp_get_kind(const int nodes, const int degree, const int* num_degrees, const int symmetries,
			 const int procs, const bool is_cpu);
extern double apsp_get_mem_usage(const int kind, const int nodes, const int degree, const int symmetries,
                                 const int *num_degrees, const int procs, const bool is_cpu);
extern void apsp_matmul(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
			const int *restrict num_degrees, const int *restrict adjacency, const bool enable_avx2, const int elements, const int symmetries);
extern void apsp_matmul_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
			      const int *restrict num_degrees, const int *restrict adjacency, const bool enable_avx2, const int symmetries);
extern void apsp_malloc(uint64_t **a, const size_t s, const bool enable_avx2);
extern void apsp_free(uint64_t *a, const bool enable_avx2);

static void apsp_mat(const int* restrict adjacency,
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
    apsp_matmul(_A, _B, _nodes, _degree, _num_degrees, adjacency, 
		_enable_avx2, _elements, _symmetries);
    
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

static void apsp_mat_saving(const int* restrict adjacency,
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
      apsp_matmul_CHUNK(_A, _B, _nodes, _degree, _num_degrees, adjacency, _enable_avx2, _symmetries);

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

void apsp_init_s(const int nodes, const int degree, const int* num_degrees, const int symmetries)
{
  _start_t = time(NULL);
  
  if(nodes % symmetries != 0)
    ERROR("nodes(%d) must be divisible by symmetries(%d)\n", nodes, symmetries);
  else if(CPU_CHUNK % 4 != 0)
    ERROR("CPU_CHUNK(%d) in parameter.h must be multiple of 4\n", CPU_CHUNK);

  _kind = apsp_get_kind(nodes, degree, num_degrees, symmetries, 1, true);
  _mem_usage = apsp_get_mem_usage(_kind, nodes, degree, symmetries, num_degrees, 1, true);
  _elements = (nodes/symmetries+(UINT64_BITS-1))/UINT64_BITS;
#ifdef __AVX2__
  if(_elements >= 4){ // For performance
    _enable_avx2 = true;
    _elements = ((_elements+3)/4)*4;  // _elements must be multiple of 4
  }
#endif

  size_t s = (_kind == APSP_NORMAL)? _elements : CPU_CHUNK;
  apsp_malloc(&_A, nodes*s*sizeof(uint64_t), _enable_avx2); // uint64_t A[nodes][s];
  apsp_malloc(&_B, nodes*s*sizeof(uint64_t), _enable_avx2); // uint64_t B[nodes][s];

  _nodes = nodes;
  _degree = degree;
  _num_degrees = num_degrees;
  _symmetries = symmetries;
  _is_profile = apsp_check_profile();
  _elapsed_time = 0;
  _times = 0;
}

void apsp_init(const int nodes, const int degree, const int* num_degrees)
{
  apsp_init_s(nodes, degree, num_degrees, 1);
}

void apsp_finalize()
{
  apsp_free(_A, _enable_avx2);
  apsp_free(_B, _enable_avx2);

  if(_is_profile){
#ifdef _OPENMP
    apsp_profile("THREADS", _kind, _symmetries, _mem_usage,
		 _start_t, time(NULL), _elapsed_time, _times, 1);
#else
    apsp_profile("SERIAL",  _kind, _symmetries, _mem_usage,
		 _start_t, time(NULL), _elapsed_time, _times, 1);
#endif
  }
}

void apsp_run(const int* restrict adjacency, int *diameter, long *sum, double *ASPL)
{
  double t = apsp_get_time();
  
  if(_kind == APSP_NORMAL)
    apsp_mat       (adjacency, diameter, sum, ASPL);
  else
    apsp_mat_saving(adjacency, diameter, sum, ASPL);

  _elapsed_time += apsp_get_time() - t;
  
  if(*diameter > _nodes)
    ERROR("This graph is not connected graph.\n");
  
  _times++;
}
