#include "common.h"

#ifdef __AVX2__
static void matmul_nregular_avx2(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                                 const int *restrict num_degrees, const int *restrict adjacency, const int elements,
                                 const int quarter_elements)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    __m256i *b = (__m256i *)(B + i*elements);
    for(int j=0;j<num_degrees[i];j++){
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

static void matmul_regular_avx2(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                                const int *restrict adjacency, const int elements, const int quarter_elements)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    __m256i *b = (__m256i *)(B + i*elements);
    for(int j=0;j<degree;j++){
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
#endif

static void matmul_nregular(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                            const int *num_degrees, const int *restrict adjacency,  const int elements)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    for(int j=0;j<num_degrees[i];j++){
      int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
      for(int k=0;k<elements;k++)
        B[i*elements+k] |= A[n*elements+k];
    }
  }
}

static void matmul_regular(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                           const int *restrict adjacency,  const int elements)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    for(int j=0;j<degree;j++){
      int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
      for(int k=0;k<elements;k++)
        B[i*elements+k] |= A[n*elements+k];
    }
  }
}

#ifdef __AVX2__
static void matmul_nregular_avx2_s(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                                   const int *restrict num_degrees, const int *restrict adjacency, const int elements,
                                   const int quarter_elements, const int based_nodes)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    __m256i *b = (__m256i *)(B + i*elements);
    int p = i/based_nodes;
    int m = i - p * based_nodes;
    for(int j=0;j<num_degrees[i];j++){
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

static void matmul_regular_avx2_s(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                                  const int *restrict adjacency, const int elements, const int quarter_elements, const int based_nodes)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    __m256i *b = (__m256i *)(B + i*elements);
    int p = i/based_nodes;
    int m = i - p * based_nodes;
    for(int j=0;j<degree;j++){
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
#endif

static void matmul_nregular_s(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                              const int *num_degrees, const int *restrict adjacency,  const int elements, const int based_nodes)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    int p = i/based_nodes;
    int m = i - p * based_nodes;
    for(int j=0;j<num_degrees[i];j++){
      int n = *(adjacency + m * degree + j) + p * based_nodes;
      if(n >= nodes) n -= nodes;
      for(int k=0;k<elements;k++)
        B[i*elements+k] |= A[n*elements+k];
    }
  }
}

static void matmul_regular_s(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                             const int *restrict adjacency,  const int elements, const int based_nodes)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    int p = i/based_nodes;
    int m = i - p * based_nodes;
    for(int j=0;j<degree;j++){
      int n = *(adjacency + m * degree + j) + p * based_nodes;
      if(n >= nodes) n -= nodes;
      for(int k=0;k<elements;k++)
        B[i*elements+k] |= A[n*elements+k];
    }
  }
}

void ODP_Matmul(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                const int *restrict num_degrees, const int *restrict adjacency, const bool enable_avx2,
                const int elements, const int symmetries)
{
#ifdef __AVX2__
  if(symmetries == 1){
    if(enable_avx2){
      if(num_degrees) matmul_nregular_avx2(A, B, nodes, degree, num_degrees, adjacency, elements, elements/4);
      else            matmul_regular_avx2(A, B, nodes, degree, adjacency, elements, elements/4);
    }
    else{
      if(num_degrees) matmul_nregular(A, B, nodes, degree, num_degrees, adjacency, elements);
      else            matmul_regular(A, B, nodes, degree, adjacency, elements);
    }
  }
  else{
    if(enable_avx2){
      if(num_degrees) matmul_nregular_avx2_s(A, B, nodes, degree, num_degrees, adjacency, elements, elements/4, nodes/symmetries);
      else            matmul_regular_avx2_s(A, B, nodes, degree, adjacency, elements, elements/4, nodes/symmetries);
    }
    else{
      if(num_degrees) matmul_nregular_s(A, B, nodes, degree, num_degrees, adjacency, elements, nodes/symmetries);
      else            matmul_regular_s(A, B, nodes, degree, adjacency, elements, nodes/symmetries);
    }
  }
#else
  if(symmetries == 1){
    if(num_degrees) matmul_nregular(A, B, nodes, degree, num_degrees, adjacency, elements);
    else            matmul_regular(A, B, nodes, degree, adjacency, elements);
  }
  else{
    if(num_degrees) matmul_nregular_s(A, B, nodes, degree, num_degrees, adjacency, elements, nodes/symmetries);
    else            matmul_regular_s(A, B, nodes, degree, adjacency, elements, nodes/symmetries);
  }
#endif
}

#ifdef __AVX2__
static void matmul_nregular_avx2_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                                       const int *num_degrees, const int *restrict adjacency)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    __m256i *b = (__m256i *)(B + i*CPU_CHUNK);
    for(int j=0;j<num_degrees[i];j++){
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

static void matmul_regular_avx2_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                                      const int *restrict adjacency)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    __m256i *b = (__m256i *)(B + i*CPU_CHUNK);
    for(int j=0;j<degree;j++){
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
#endif

static void matmul_nregular_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                                  const int *num_degrees, const int *restrict adjacency)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    for(int j=0;j<num_degrees[i];j++){
      int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
      for(int k=0;k<CPU_CHUNK;k++)
        B[i*CPU_CHUNK+k] |= A[n*CPU_CHUNK+k];
    }
  }
}

static void matmul_regular_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                                 const int *restrict adjacency)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    for(int j=0;j<degree;j++){
      int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
      for(int k=0;k<CPU_CHUNK;k++)
        B[i*CPU_CHUNK+k] |= A[n*CPU_CHUNK+k];
    }
  }
}

#ifdef __AVX2__
static void matmul_nregular_avx2_CHUNK_s(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                                         const int *num_degrees, const int *restrict adjacency, const int based_nodes)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    __m256i *b = (__m256i *)(B + i*CPU_CHUNK);
    int p = i/based_nodes;
    int m = i - p * based_nodes;
    for(int j=0;j<num_degrees[i];j++){
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

static void matmul_regular_avx2_CHUNK_s(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                                        const int *restrict adjacency, const int based_nodes)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    __m256i *b = (__m256i *)(B + i*CPU_CHUNK);
    int p = i/based_nodes;
    int m = i - p * based_nodes;
    for(int j=0;j<degree;j++){
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
#endif

static void matmul_nregular_CHUNK_s(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                                    const int *num_degrees, const int *restrict adjacency, const int based_nodes)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    int p = i/based_nodes;
    int m = i - p * based_nodes;
    for(int j=0;j<num_degrees[i];j++){
      int n = *(adjacency + m * degree + j) + p * based_nodes;
      if(n >= nodes) n -= nodes;
      for(int k=0;k<CPU_CHUNK;k++)
        B[i*CPU_CHUNK+k] |= A[n*CPU_CHUNK+k];
    }
  }
}

static void matmul_regular_CHUNK_s(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                                   const int *restrict adjacency, const int based_nodes)
{
#pragma omp parallel for
  for(int i=0;i<nodes;i++){
    int p = i/based_nodes;
    int m = i - p * based_nodes;
    for(int j=0;j<degree;j++){
      int n = *(adjacency + m * degree + j) + p * based_nodes;
      if(n >= nodes) n -= nodes;
      for(int k=0;k<CPU_CHUNK;k++)
        B[i*CPU_CHUNK+k] |= A[n*CPU_CHUNK+k];
    }
  }
}

void ODP_Matmul_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
                      const int *num_degrees, const int *restrict adjacency, const bool enable_avx2, const int symmetries)
{
#ifdef __AVX2__
  if(symmetries == 1){
    if(enable_avx2){
      if(num_degrees) matmul_nregular_avx2_CHUNK(A, B, nodes, degree, num_degrees, adjacency);
      else            matmul_regular_avx2_CHUNK(A, B, nodes, degree, adjacency);
    }
    else{
      if(num_degrees) matmul_nregular_CHUNK(A, B, nodes, degree, num_degrees, adjacency);
      else            matmul_regular_CHUNK(A, B, nodes, degree, adjacency);
    }
  }
  else{ // symmetries != 1
    if(enable_avx2){
      if(num_degrees) matmul_nregular_avx2_CHUNK_s(A, B, nodes, degree, num_degrees, adjacency, nodes/symmetries);
      else            matmul_regular_avx2_CHUNK_s(A, B, nodes, degree, adjacency, nodes/symmetries);
    }
    else{
      if(num_degrees) matmul_nregular_CHUNK_s(A, B, nodes, degree, num_degrees, adjacency, nodes/symmetries);
      else            matmul_regular_CHUNK_s(A, B, nodes, degree, adjacency, nodes/symmetries);
    }
  }
#endif
  if(symmetries == 1){
    if(num_degrees) matmul_nregular_CHUNK(A, B, nodes, degree, num_degrees, adjacency);
    else            matmul_regular_CHUNK(A, B, nodes, degree, adjacency);
  }
  else{ // symmetries != 1
    if(num_degrees) matmul_nregular_CHUNK_s(A, B, nodes, degree, num_degrees, adjacency, nodes/symmetries);
    else            matmul_regular_CHUNK_s(A, B, nodes, degree, adjacency, nodes/symmetries);
  }
}
