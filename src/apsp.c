#include "common.h"
static uint64_t *_A, *_B;
static int _nodes, _degree, _groups, _kind;
static const int* _num_degrees;
static double _mem_usage;

extern void apsp_start_profile();
extern void apsp_end_profile(const char* name, const int kind, const int groups, const double mem_usage, const int procs);
extern int  apsp_get_kind(const int nodes, const int degree, const int* num_degrees, const int groups,
			  const int procs);
extern double apsp_get_mem_usage(const int kind, const int nodes, const int degree, const int groups,
                                 const int *num_degree, const int procs);

static void apsp_mat(const int* restrict adjacency,
		     int *diameter, long *sum, double *ASPL)
{
  unsigned int elements = (_nodes/_groups+(UINT64_BITS-1))/UINT64_BITS;

#pragma omp parallel for
  for(int i=0;i<_nodes*elements;i++)
    _A[i] = _B[i] = 0;

#pragma omp parallel for
  for(int i=0;i<_nodes/_groups;i++){
    unsigned int offset = i*elements+i/UINT64_BITS;
    _A[offset] = _B[offset] = (0x1ULL << (i%UINT64_BITS));
  }

  *sum = (long)_nodes * (_nodes - 1);
  *diameter = 1;
  for(int kk=0;kk<_nodes;kk++){
    if(!_num_degrees){
#pragma omp parallel for
      for(int i=0;i<_nodes;i++)
	for(int j=0;j<_degree;j++){
	  int n = *(adjacency + i * _degree + j);  // int n = adjacency[i][j];
	  for(int k=0;k<elements;k++)
	    _B[i*elements+k] |= _A[n*elements+k];
	}
    }
    else{
#pragma omp parallel for
      for(int i=0;i<_nodes;i++)
        for(int j=0;j<_num_degrees[i];j++){
          int n = *(adjacency + i * _degree + j);  // int n = adjacency[i][j];
          for(int k=0;k<elements;k++)
            _B[i*elements+k] |= _A[n*elements+k];
        }
    }

    uint64_t num = 0;
#pragma omp parallel for reduction(+:num)
    for(int i=0;i<elements*_nodes;i++)
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
  unsigned int elements = (_nodes/_groups+(UINT64_BITS-1))/UINT64_BITS;
  int parsize = (elements+(CHUNK-1))/CHUNK;

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

  *ASPL = *sum / (((double)_nodes-1)*_nodes);
  *sum /= 2.0;
}

void apsp_init_s(const int nodes, const int degree, const int* num_degrees, const int groups)
{
  if(nodes % groups != 0)
    ERROR("nodes(%d) must be divisible by group(%d)\n", nodes, groups);

  _kind = apsp_get_kind(nodes, degree, num_degrees, groups, 1);
  _mem_usage = apsp_get_mem_usage(_kind, nodes, degree, groups, num_degrees, 1);
  if(_kind == APSP_NORMAL){
    unsigned int elements = (nodes/groups+(UINT64_BITS-1))/UINT64_BITS;
    size_t s = nodes * elements * sizeof(uint64_t);
    posix_memalign((void **)&_A, ALIGN_VALUE, s); // uint64_t A[nodes][elements];
    posix_memalign((void **)&_B, ALIGN_VALUE, s); // uint64_t B[nodes][elements];
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
}

void apsp_init(const int nodes, const int degree, const int* num_degrees)
{
  apsp_init_s(nodes, degree, num_degrees, 1);
}

void apsp_finalize()
{
  free(_A);
  free(_B);
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
