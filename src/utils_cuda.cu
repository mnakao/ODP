#include "common.h"

__global__ void create_adjacency(const int based_elements, const int total_elements,
				 const int based_nodes, const int nodes, int* __restrict__ adjacency_dev)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x + based_elements;
  while (tid < total_elements) {
    int t = tid/based_elements;
    int i = tid - (t*based_elements);
    int v = adjacency_dev[i] + t*based_nodes;
    adjacency_dev[tid] = (v < nodes)? v : v - nodes;

    tid += blockDim.x * gridDim.x;
  }
}

__global__ void clear_buffers(uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int length)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid<length) {
    A[tid] = B[tid] = 0;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void popcnt(const uint64_t* __restrict__ B, const int nodes,
		       const unsigned int elements, uint64_t* __restrict__ result)
{
  __shared__ uint64_t cache[THREADS];
  int cacheIndex = threadIdx.x;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  uint64_t num = 0;
  while (tid < elements*nodes) {
    num += POPCNT(B[tid]);
    tid += blockDim.x * gridDim.x;
  }
  cache[cacheIndex] = num;
  __syncthreads();

  int i = blockDim.x/2;
  while (i != 0){
    if (cacheIndex < i)
      cache[cacheIndex] += cache[cacheIndex+i];
    __syncthreads();
    i /= 2;
  }

  if(cacheIndex == 0)
    result[blockIdx.x] = cache[0];
}

__global__ void matrix_op(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
			  const int* __restrict__ num_degrees, const int nodes, const int degree, const unsigned int elements)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if(!num_degrees){
    while (tid < nodes*elements) {
      int i = tid / elements;
      int k = tid % elements;
      uint64_t tmp = B[tid];
      for(int j=0;j<degree;j++){
        int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
        tmp |= A[n*elements+k];
      }
      B[tid] = tmp;
      tid += blockDim.x * gridDim.x;
    }
  }
  else{
    while (tid < nodes*elements) {
      int i = tid / elements;
      int k = tid % elements;
      uint64_t tmp = B[tid];
      for(int j=0;j<num_degrees[i];j++){
        int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
        tmp |= A[n*elements+k];
      }
      B[tid] = tmp;
      tid += blockDim.x * gridDim.x;
    }
  }
}

