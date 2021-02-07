#include "common.h"
extern int ODP_GLOBAL_ADJ_GRID(const int width, const int height, const int degree, const int symmetries,
			       const int *adjacency, const int v, const int d);

static void matmul(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
		   const int *num_degrees, const int *restrict adjacency,  const int elements)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    int d = (!num_degrees)? degree : num_degrees[i];
    for(int j=0;j<d;j++){
      int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
      for(int k=0;k<elements;k++)
        B[i*elements+k] |= A[n*elements+k];
    }
  }
}

static void matmul_s(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int height,
                     const int degree, const int *num_degrees, const int *restrict adjacency, const int elements,
                     const int symmetries, const int itable[nodes])
{
  int based_nodes = nodes/symmetries;
  if(!itable){
#pragma omp parallel for
    for(int i=0;i<nodes;i++){
      int p = i/based_nodes;
      int m = i - p * based_nodes;
      int d = (!num_degrees)? degree : num_degrees[i];
      for(int j=0;j<d;j++){
        int n = *(adjacency + m * degree + j) + p * based_nodes;
        if(n >= nodes) n -= nodes;
        for(int k=0;k<elements;k++)
          B[i*elements+k] |= A[n*elements+k];
      }
    }
  }
  else{
    int width = nodes/height;
#pragma omp parallel for
    for(int i=0;i<nodes;i++){
      int ii = itable[i];
      int d = (!num_degrees)? degree : num_degrees[i];
      for(int j=0;j<d;j++){
        int n = ODP_GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, i, j);
        int nn = itable[n];
        for(int k=0;k<elements;k++)
          B[ii*elements+k] |= A[nn*elements+k];
      }
    }
  }
}

#ifdef __AVX2__
static void matmul_avx2(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                        const int *restrict num_degrees, const int *restrict adjacency, const int elements,
                        const int quarter_elements)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    __m256i *b = (__m256i *)(B + i*elements);
    int d = (!num_degrees)? degree : num_degrees[i];
    for(int j=0;j<d;j++){
      int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
      __m256i *a = (__m256i *)(A + n*elements);
      for(int k=0;k<quarter_elements;k++){
        __m256i aa = _mm256_load_si256(a+k);
        __m256i bb = _mm256_load_si256(b+k);
        _mm256_store_si256(b+k, _mm256_or_si256(aa, bb));
      }
    }
  }
}

static void matmul_avx2_s(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int height, const int degree,
			  const int *restrict num_degrees, const int *restrict adjacency, const int elements,
			  const int quarter_elements, const int symmetries, const int itable[nodes])
{
  int based_nodes = nodes/symmetries;
  if(!itable){
#pragma omp parallel for
    for(int i=0;i<nodes;i++){
      __m256i *b = (__m256i *)(B + i*elements);
      int p = i/based_nodes;
      int m = i - p * based_nodes;
      int d = (!num_degrees)? degree : num_degrees[i];
      for(int j=0;j<d;j++){
	int n = *(adjacency + m * degree + j) + p * based_nodes;
	if(n >= nodes) n -= nodes;
	__m256i *a = (__m256i *)(A + n*elements);
	for(int k=0;k<quarter_elements;k++){
	  __m256i aa = _mm256_load_si256(a+k);
	  __m256i bb = _mm256_load_si256(b+k);
	  _mm256_store_si256(b+k, _mm256_or_si256(aa, bb));
	}
      }
    }
  }
  else{
    int width = nodes/height;
#pragma omp parallel for
    for(int i=0;i<nodes;i++){
      __m256i *b = (__m256i *)(B + itable[i]*elements);
      int d = (!num_degrees)? degree : num_degrees[i];
      for(int j=0;j<d;j++){
	int n = ODP_GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, i, j);
	__m256i *a = (__m256i *)(A + itable[n]*elements);
	for(int k=0;k<quarter_elements;k++){
	  __m256i aa = _mm256_load_si256(a+k);
	  __m256i bb = _mm256_load_si256(b+k);
	  _mm256_store_si256(b+k, _mm256_or_si256(aa, bb));
	}
      }
    }
  }
}
#endif

void ODP_Matmul(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int height, const int degree,
                const int *restrict num_degrees, const int *restrict adjacency, const bool enable_avx2,
                const int elements, const int symmetries, const int itable[nodes])
{
#ifdef __AVX2__
  if(symmetries == 1){
    if(enable_avx2){
      matmul_avx2(A, B, nodes, degree, num_degrees, adjacency, elements, elements/4);
    }
    else{
      matmul(A, B, nodes, degree, num_degrees, adjacency, elements);
    }
  }
  else{
    if(enable_avx2){
      matmul_avx2_s(A, B, nodes, height, degree, num_degrees, adjacency, elements, elements/4, symmetries, itable);
    }
    else{
      matmul_s(A, B, nodes, height, degree, num_degrees, adjacency, elements, symmetries, itable);
    }
  }
#else
  if(symmetries == 1){
    matmul(A, B, nodes, degree, num_degrees, adjacency, elements);
  }
  else{
    matmul_s(A, B, nodes, height, degree, num_degrees, adjacency, elements, symmetries, itable);
  }
#endif
}

static void matmul_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                         const int *num_degrees, const int *restrict adjacency)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    int d = (!num_degrees)? degree : num_degrees[i];
    for(int j=0;j<d;j++){
      int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
      for(int k=0;k<CPU_CHUNK;k++)
        B[i*CPU_CHUNK+k] |= A[n*CPU_CHUNK+k];
    }
  }
}

static void matmul_CHUNK_s(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int height,
                           const int degree, const int *num_degrees, const int *restrict adjacency,
                           const int symmetries, const int itable[nodes])
{
  int based_nodes = nodes/symmetries;
  if(!itable){
#pragma omp parallel for
    for(int i=0;i<nodes;i++){
      int p = i/based_nodes;
      int m = i - p * based_nodes;
      int d = (!num_degrees)? degree : num_degrees[i];
      for(int j=0;j<d;j++){
        int n = *(adjacency + m * degree + j) + p * based_nodes;
        if(n >= nodes) n -= nodes;
        for(int k=0;k<CPU_CHUNK;k++)
          B[i*CPU_CHUNK+k] |= A[n*CPU_CHUNK+k];
      }
    }
  }
  else{
    int width = nodes/height;
#pragma omp parallel for
    for(int i=0;i<nodes;i++){
      int ii = itable[i];
      int d = (!num_degrees)? degree : num_degrees[i];
      for(int j=0;j<d;j++){
        int n = ODP_GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, i, j);
        int nn = itable[n];
        for(int k=0;k<CPU_CHUNK;k++)
          B[ii*CPU_CHUNK+k] |= A[nn*CPU_CHUNK+k];
      }
    }
  }
}

#ifdef __AVX2__
static void matmul_avx2_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
			      const int *num_degrees, const int *restrict adjacency)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    __m256i *b = (__m256i *)(B + i*CPU_CHUNK);
    int d = (!num_degrees)? degree : num_degrees[i];
    for(int j=0;j<d;j++){
      int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
      __m256i *a = (__m256i *)(A + n*CPU_CHUNK);
      for(int k=0;k<CPU_CHUNK/4;k++){
        __m256i aa = _mm256_load_si256(a+k);
        __m256i bb = _mm256_load_si256(b+k);
        _mm256_store_si256(b+k, _mm256_or_si256(aa, bb));
      }
    }
  }
}

static void matmul_avx2_CHUNK_s(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int height,
				const int degree, const int *num_degrees, const int *restrict adjacency, const int symmetries,
				const int itable[nodes])
{
  int based_nodes = nodes/symmetries;
  if(!itable){
#pragma omp parallel for
    for(int i=0;i<nodes;i++){
      __m256i *b = (__m256i *)(B + i*CPU_CHUNK);
      int p = i/based_nodes;
      int m = i - p * based_nodes;
      int d = (!num_degrees)? degree : num_degrees[i];
      for(int j=0;j<d;j++){
	int n = *(adjacency + m * degree + j) + p * based_nodes;
	if(n >= nodes) n -= nodes;
	__m256i *a = (__m256i *)(A + n*CPU_CHUNK);
	for(int k=0;k<CPU_CHUNK/4;k++){
	  __m256i aa = _mm256_load_si256(a+k);
	  __m256i bb = _mm256_load_si256(b+k);
	  _mm256_store_si256(b+k, _mm256_or_si256(aa, bb));
	}
      }
    }
  }
  else{
    int width = nodes/height;
#pragma omp parallel for
    for(int i=0;i<nodes;i++){
      __m256i *b = (__m256i *)(B + itable[i]*CPU_CHUNK);
      int d = (!num_degrees)? degree : num_degrees[i];
      for(int j=0;j<d;j++){
	int n = ODP_GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, i, j);
        __m256i *a = (__m256i *)(A + itable[n]*CPU_CHUNK);
        for(int k=0;k<CPU_CHUNK/4;k++){
          __m256i aa = _mm256_load_si256(a+k);
          __m256i bb = _mm256_load_si256(b+k);
          _mm256_store_si256(b+k, _mm256_or_si256(aa, bb));
        }
      }
    }
  }
}
#endif

void ODP_Matmul_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int height, const int degree,
                      const int *num_degrees, const int *restrict adjacency, const bool enable_avx2, const int symmetries,
		      const int itable[nodes])
{
#ifdef __AVX2__
  if(symmetries == 1){
    if(enable_avx2){
      matmul_avx2_CHUNK(A, B, nodes, degree, num_degrees, adjacency);
    }
    else{
      matmul_CHUNK(A, B, nodes, degree, num_degrees, adjacency);
    }
  }
  else{ // symmetries != 1
    if(enable_avx2){
      matmul_avx2_CHUNK_s(A, B, height, nodes, degree, num_degrees, adjacency, symmetries, itable);
    }
    else{
      matmul_CHUNK_s(A, B, nodes, height, degree, num_degrees, adjacency, symmetries, itable);
    }
  }
#else
  if(symmetries == 1){
    matmul_CHUNK(A, B, nodes, degree, num_degrees, adjacency);
  }
  else{ // symmetries != 1
    matmul_CHUNK_s(A, B, nodes, height, degree, num_degrees, adjacency, symmetries, itable);
  }
#endif
}
