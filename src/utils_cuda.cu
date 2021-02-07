#include "common.h"
#define WIDTH(v,h) ((v)/(h))
#define HEIGHT(v,h) ((v)%(h))

__device__ static int rotate(const int v, const int width, const int height, const int symmetries, const int degree)
{
  int w = WIDTH (v, height);
  int h = HEIGHT(v, height);
  if(symmetries == 2)
    return (width-w-1)*height + (height-h-1);

  if(degree == 90)       return h*height + (height-w-1);
  else if(degree == 180) return (height-w-1)*height + (height-h-1);
  else                   return (height-h-1)*height + w; // degree == 270
}

__device__ static int local_index(const int x, const int width, const int height, const int symmetries)
{
  if(symmetries == 2){
    return x;
  }
  else{ // symmetries == 4
    int based_height = height/2;
    return WIDTH(x,height)*based_height + HEIGHT(x,height);
  }
}

__device__ static int global_adj(const int width, const int height, const int degree, const int symmetries,
				 const int *adjacency, const int v, const int d)
{
  if(symmetries == 1){
    return adjacency[v*degree+d];
  }
  else if(symmetries == 2){
    int based_width = width/2;
    if(WIDTH(v,height) < based_width)
      return adjacency[v*degree+d];
    else{
      int y = adjacency[rotate(v, width, height, symmetries, 180)*degree+d];
      return rotate(y, width, height, symmetries, 180);
    }
  }
  else{ // symmetries == 4
    int based_width  = width/2;
    int based_height = height/2;
    if(WIDTH(v,height) < based_width && HEIGHT(v,height) < based_height){
      return adjacency[local_index(v,width,height,symmetries)*degree+d];
    }
    else if(WIDTH(v,height) < based_width && HEIGHT(v,height) >= based_height){
      int x = rotate(v, width, height, symmetries, 270);
      int y = adjacency[local_index(x,width,height,symmetries)*degree+d];
      return rotate(y, width, height, symmetries, 90);
    }
    else if(WIDTH(v,height) >= based_width && HEIGHT(v,height) >= based_height){
      int x = rotate(v, width, height, symmetries, 180);
      int y = adjacency[local_index(x,width,height,symmetries)*degree+d];
      return rotate(y, width, height, symmetries, 180);
    }
    else{
      int x = rotate(v, width, height, symmetries, 90);
      int y = adjacency[local_index(x,width,height,symmetries)*degree+d];
      return rotate(y, width, height, symmetries, 270);
    }
  }
}

__global__ void ODP_Clear_buffers(uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int length)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid<length) {
    A[tid] = B[tid] = 0;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void ODP_Popcnt(const uint64_t* __restrict__ B, const int nodes,
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

__global__ void ODP_Matmul_cuda(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
				const int* __restrict__ num_degrees, const int nodes, const int height, const int degree,
				const unsigned int elements, const int symmetries, const int* __restrict__ itable)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(symmetries == 1){
    while (tid < nodes*elements) {
      int i = tid / elements;
      int k = tid % elements;
      uint64_t tmp = B[tid];
      int d = (!num_degrees)? degree : num_degrees[i];
      for(int j=0;j<d;j++){
	int n = *(adjacency + i * degree + j);
	tmp |= A[n*elements+k];
      }
      B[tid] = tmp;
      tid += blockDim.x * gridDim.x;
    }
  }
  else{
    int based_nodes = nodes/symmetries;
    if(!itable){
      while (tid < nodes*elements) {
	int i = tid / elements;
	int k = tid % elements;
	int t = i / based_nodes;
	int h = t * based_nodes;
	int m = i - h;
	uint64_t tmp = B[tid];
	int d = (!num_degrees)? degree : num_degrees[i];
	for(int j=0;j<d;j++){
	  int n = *(adjacency + m * degree + j) + h; // int n = adjacency[i][j] <- adjacency[m][j] + t*based_nodes;
	  if(n >= nodes) n -= nodes;
	  tmp |= A[n*elements+k];
	}
	B[tid] = tmp;
	tid += blockDim.x * gridDim.x;
      }
    }
    else{
      int width = nodes/height;
      while (tid < nodes*elements) {
	int i = tid / elements;
	int k = tid % elements;
	int ii = itable[i]*elements+k;
	uint64_t tmp = B[ii];
	int d = (!num_degrees)? degree : num_degrees[i];
	for(int j=0;j<d;j++){
	  int n = global_adj(width, height, degree, symmetries, adjacency, i, j);
	  tmp |= A[itable[n]*elements+k];
	}
	B[ii] = tmp;
	tid += blockDim.x * gridDim.x;
      }
    }
  }
}

__global__ void ODP_Matmul_CHUNK_cuda(const uint64_t* __restrict__ A, uint64_t* __restrict__ B, const int* __restrict__ adjacency,
				      const int* __restrict__ num_degrees, const int nodes, const int height, const int degree,
				      const int symmetries, const int* __restrict__ itable)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(symmetries == 1){
    while (tid < nodes*GPU_CHUNK) {
      int i = tid / GPU_CHUNK;
      int k = tid % GPU_CHUNK;
      uint64_t tmp = B[tid];
      int d = (!num_degrees)? degree : num_degrees[i];
      for(int j=0;j<d;j++){
        int n = *(adjacency + i * degree + j);
        tmp |= A[n*GPU_CHUNK+k];
      }
      B[tid] = tmp;
      tid += blockDim.x * gridDim.x;
    }
  }
  else{
    int based_nodes = nodes/symmetries;
    if(!itable){
      while (tid < nodes*GPU_CHUNK) {
	int i = tid / GPU_CHUNK;
	int k = tid % GPU_CHUNK;
	int t = i / based_nodes;
	int h = t * based_nodes;
	int m = i - h;
	uint64_t tmp = B[tid];
	int d = (!num_degrees)? degree : num_degrees[i];
	for(int j=0;j<d;j++){
	  int n = *(adjacency + m * degree + j) + h; // int n = adjacency[i][j] <- adjacency[m][j] + t*based_nodes;
	  if(n >= nodes) n -= nodes;
	  tmp |= A[n*GPU_CHUNK+k];
	}
	B[tid] = tmp;
	tid += blockDim.x * gridDim.x;
      }
    }
    else{
      int width = nodes/height;
      while (tid < nodes*GPU_CHUNK) {
	int i = tid / GPU_CHUNK;
	int k = tid % GPU_CHUNK;
	int ii = itable[i]*GPU_CHUNK+k;
	uint64_t tmp = B[ii];
	int d = (!num_degrees)? degree : num_degrees[i];
	for(int j=0;j<d;j++){
	  int n = global_adj(width, height, degree, symmetries, adjacency, i, j);
	  tmp |= A[itable[n]*GPU_CHUNK+k];
	}
	B[ii] = tmp;
	tid += blockDim.x * gridDim.x;
      }
    }
  }
}
