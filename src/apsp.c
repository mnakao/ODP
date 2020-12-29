#include "common.h"
static uint64_t *_A, *_B;
static int _nodes, _degree, _groups, _kind;
static const int* _num_degrees;
static unsigned int _elements;
static double _mem_usage;

extern void apsp_start_profile();
extern void apsp_end_profile(const char* name, const int kind, const int groups, const double mem_usage, const int procs);
extern int  apsp_get_kind(const int nodes, const int degree, const int* num_degrees, const int groups,
			  const int procs);
extern double apsp_get_mem_usage(const int kind, const int nodes, const int degree, const int groups,
                                 const int *num_degree, const int procs);
extern void matmul(const uint64_t *restrict A, uint64_t *restrict B, const int* restrict adjacency,
		   const int nodes, int degree, const int *num_degrees, const unsigned int rows);

static void apsp_mat(const int* restrict adjacency,
		     int *diameter, long *sum, double *ASPL)
{
#pragma omp parallel for
  for(int i=0;i<_nodes*_elements;i++)
    _A[i] = _B[i] = 0;

#pragma omp parallel for
  for(int i=0;i<_nodes/_groups;i++){
    unsigned int offset = i*_elements+i/UINT64_BITS;
    _A[offset] = _B[offset] = (0x1ULL << (i%UINT64_BITS));
  }

  *sum = (long)_nodes * (_nodes - 1);
  *diameter = 1;
  for(int kk=0;kk<_nodes;kk++){
    matmul(_A, _B, adjacency, _nodes, _degree, _num_degrees, _elements);

    uint64_t num = 0;
#pragma omp parallel for reduction(+:num)
    for(int i=0;i<_elements*_nodes;i++)
      num += POPCNT(_B[i]);

    num *= _groups;
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
  int parsize = (_elements+(CHUNK-1))/CHUNK;

  *sum = (long)_nodes * (_nodes - 1);
  *diameter = 1;
  for(int t=0;t<parsize;t++){
    unsigned int kk, l;
#pragma omp parallel for
    for(int i=0;i<_nodes*CHUNK;i++)
      _A[i] = _B[i] = 0;
    
    for(l=0; l<UINT64_BITS*CHUNK && UINT64_BITS*t*CHUNK+l<_nodes/_groups; l++){
      unsigned int offset = (UINT64_BITS*t*CHUNK+l)*CHUNK+l/UINT64_BITS;
      _A[offset] = _B[offset] = (0x1ULL<<(l%UINT64_BITS));
    }

    for(kk=0;kk<_nodes;kk++){
      matmul(_A, _B, adjacency, _nodes, _degree, _num_degrees, CHUNK);

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

  *ASPL = *sum / (((double)_nodes-1)*_nodes);
  *sum /= 2.0;
}

void apsp_init_s(const int nodes, const int degree, const int* num_degrees, const int groups)
{
  if(nodes % groups != 0)
    ERROR("nodes(%d) must be divisible by group(%d)\n", nodes, groups);
  else if(CHUNK % 4 != 0)
    ERROR("CHUNK(%d) in parameter.h must be multiple of 4\n", CHUNK);

  _kind = apsp_get_kind(nodes, degree, num_degrees, groups, 1);
  _mem_usage = apsp_get_mem_usage(_kind, nodes, degree, groups, num_degrees, 1);
  _elements = (nodes/groups+(UINT64_BITS-1))/UINT64_BITS;
#ifdef __AVX2__
  _elements = ((_elements+3)/4)*4;  // _elements must be multiple of 4
#endif
  
  if(_kind == APSP_NORMAL){
    size_t s = nodes * _elements * sizeof(uint64_t);
#ifdef __AVX2__
    _A = (uint64_t *) _mm_malloc(s, ALIGN_VALUE); // uint64_t A[nodes][elements];
    _B = (uint64_t *) _mm_malloc(s, ALIGN_VALUE); // uint64_t B[nodes][elements];
#else
    posix_memalign((void **)&_A, ALIGN_VALUE, s); // uint64_t A[nodes][elements];
    posix_memalign((void **)&_B, ALIGN_VALUE, s); // uint64_t B[nodes][elements];
#endif
  }
  else{
    size_t s = nodes * CHUNK * sizeof(uint64_t); // CHUNK must be multiple of 4
#ifdef __AVX2__
    _A = (uint64_t *) _mm_malloc(s, ALIGN_VALUE); // uint64_t A[nodes][CHUNK];
    _B = (uint64_t *) _mm_malloc(s, ALIGN_VALUE); // uint64_t B[nodes][CHUNK];
#else
    posix_memalign((void **)&_A, ALIGN_VALUE, s); // uint64_t A[nodes][CHUNK];
    posix_memalign((void **)&_B, ALIGN_VALUE, s); // uint64_t B[nodes][CHUNK];
#endif
  }

  _nodes = nodes;
  _degree = degree;
  _num_degrees = num_degrees;
  _groups = groups;
}

void apsp_init(const int nodes, const int degree, const int* num_degrees)
{
  apsp_init_s(nodes, degree, num_degrees, 1);
}

void apsp_finalize()
{
#ifdef __AVX2__
  _mm_free(_A);
  _mm_free(_B);
#else
  free(_A);
  free(_B);
#endif
}

void apsp_run(const int* restrict adjacency, int *diameter, long *sum, double *ASPL)
{
  apsp_start_profile();
  
  if(_kind == APSP_NORMAL)
    apsp_mat       (adjacency, diameter, sum, ASPL);
  else
    apsp_mat_saving(adjacency, diameter, sum, ASPL);
  
  if(*diameter > _nodes)
    ERROR("This graph is not connected graph.\n");

#ifdef _OPENMP
  apsp_end_profile("THREADS", _kind, _groups, _mem_usage, 1);
#else
  apsp_end_profile("SERIAL", _kind, _groups, _mem_usage, 1);
#endif
}
