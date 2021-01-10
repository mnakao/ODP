#include "common.h"
static int *_n = NULL;
static int *_d = NULL;
static int _r, _nodes, _degree;

double apsp_get_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1.0e-6 * t.tv_usec;
}

bool apsp_check_profile()
{
  char *val = getenv("APSP_PROFILE");
  if(val){
    if(atoi(val) == 1)
      return true;
  }
  return false;
}

static void check_graph_parameters(const int nodes, const int degree)
{
  if(nodes % 2 == 1 && degree % 2 == 1)
    ERROR("Nodes(%d) or Degree(%d) must be a multiple of 2.\n", nodes, degree);
}

static int get_random(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
}

static void swap(int *a, int *b)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

static bool check_multiple_vertices(const int e00, const int e01, const int e10, const int e11)
{
  return (e00 == e10 || e01 == e11 || e00 == e11 || e01 == e10);
}

void apsp_srand(const unsigned int seed)
{
  srand(seed);
}

void apsp_write_edge_general(const int lines, const int edge[lines][2], char *fname)
{
  FILE *fp = NULL;
  
  if((fp = fopen(fname, "w")) == NULL)
    ERROR("Cannot open %s\n", fname);
  
  for(int i=0;i<lines;i++)
    fprintf(fp, "%d %d\n", edge[i][0], edge[i][1]);
  
  fclose(fp);
}

void apsp_write_edge_grid(const int lines, const int height, const int edge[lines][2], char *fname)
{
  FILE *fp = NULL;

  if((fp = fopen(fname, "w")) == NULL)
    ERROR("Cannot open %s\n", fname);

  for(int i=0;i<lines;i++)
    fprintf(fp, "%d,%d %d,%d\n",
	    WIDTH(edge[i][0], height), HEIGHT(edge[i][0], height),
            WIDTH(edge[i][1], height), HEIGHT(edge[i][1], height));
    
  fclose(fp);
}

static bool check_length(const int v, const int w, const int height, const int length)
{
  int w0 = WIDTH(v,height);
  int h0 = HEIGHT(v,height);
  int w1 = WIDTH(w,height);
  int h1 = HEIGHT(w,height);
  int distance = abs(w0 - w1) + abs(h0 - h1);
  
  return (distance <= length);
}

static void simple_2opt_grid(const int height, const int length, const int lines, int edge[lines][2])
{
  while(1){
    int e0, e1;
    do{
      e0 = get_random(lines);
      e1 = get_random(lines);
    } while(check_multiple_vertices(edge[e0][0], edge[e0][1], edge[e1][0], edge[e1][1]));

    if(get_random(2) == 0){
      if(check_length(edge[e0][0], edge[e1][1], height, length) && check_length(edge[e0][1], edge[e1][0], height, length)){
	swap(&edge[e0][1], &edge[e1][1]);
	break;
      }
    }
    else{
      if(check_length(edge[e0][0], edge[e1][0], height, length) && check_length(edge[e0][1], edge[e1][1], height, length)){
	swap(&edge[e0][1], &edge[e1][0]);
	break;
      }
    }
  }
}

// Inherited from http://research.nii.ac.jp/graphgolf/c/create-lattice.c
void apsp_generate_random_grid(const int width, const int height, const int degree, const int length, int (*edge)[2])
{
  int nodes = width * height;
  check_graph_parameters(nodes, degree);

  int i = 0;
  for(int x=0;x<width/2;x++){
    for(int y=0;y<height;y++){
      for(int k=0;k<degree;k++){
        edge[i][0] = y + 2 * x * height;
        edge[i][1] = edge[i][0] + height;
        i++;
      }
    }
  }

  if(width%2 == 1){
    for(int y=0;y<height/2;y++){
      for(int k=0;k<degree;k++){
        edge[i][0] = (width - 1) * height + 2 * y;
        edge[i][1] = edge[i][0] + 1;
        i++;
      }
    }

    /* add self-loop */
    if(height%2 == 1){
      for(int k=0;k<degree/2;k++){
        edge[i][0] = edge[i][1] = nodes - 1;
        i++;
      }
    }
  }

  int lines = (nodes*degree)/2;
  for(int i=0;i<lines*GEN_GRAPH_ITERS;i++)  // Give randomness
    simple_2opt_grid(height, length, lines, edge);
}

void apsp_malloc(uint64_t **a, const size_t s, const bool enable_avx2)
{
  if(enable_avx2)
    *a = _mm_malloc(s, ALIGN_VALUE);
  else
    posix_memalign((void **)a, ALIGN_VALUE, s);
}

void apsp_free(uint64_t *a, const bool enable_avx2)
{
  if(enable_avx2)
    _mm_free(a);
  else
    free(a);
}

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

void apsp_matmul(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
		 const int *restrict num_degrees, const int *restrict adjacency, const bool enable_avx2,
		 const int elements, const int symmetries)
{
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
}

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

void apsp_matmul_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
		       const int *num_degrees, const int *restrict adjacency, const bool enable_avx2, const int symmetries)
{
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
}

double apsp_get_mem_usage(const int kind, const int nodes, const int degree, const int symmetries,
			  const int *num_degrees, const int procs, const bool is_cpu)
{
  int Mbyte = 1024*1024;
  int chunk = (is_cpu)? CPU_CHUNK : GPU_CHUNK;
  double AB_mem = (kind == APSP_NORMAL)? (nodes*((double)nodes/(4*symmetries*procs))) : (double)16*nodes*chunk;

  if(is_cpu){
    return AB_mem/Mbyte;
  }
  else{ // on GPU
    double res_mem = (double)sizeof(uint64_t)*BLOCKS;
    double adj_mem = (double)sizeof(int)*(nodes/symmetries)*degree;
    double deg_mem = (num_degrees)? (double)sizeof(int)*nodes : 0;
    return (AB_mem+res_mem+adj_mem+deg_mem)/Mbyte;
  }
}

int apsp_get_kind(const int nodes, const int degree, const int* num_degrees, const int symmetries,
		  const int procs, const int is_cpu)
{
  char *val = getenv("APSP");
  int kind;
  if(val == NULL){
    double normal_mem_usage = apsp_get_mem_usage(APSP_NORMAL, nodes, degree, symmetries, num_degrees, procs, is_cpu);
    if(normal_mem_usage <= MEM_THRESHOLD)
      kind = APSP_NORMAL;
    else
      kind = APSP_SAVING;
  }
  else if(strcmp(val, "NORMAL") == 0){
    kind = APSP_NORMAL;
  }
  else if(strcmp(val, "SAVING") == 0){
    kind = APSP_SAVING;
  }
  else{
    ERROR("Unknown APSP value (%s)\n", val);
  }

  return kind;
}

void apsp_profile(const char* name, const int kind, const int symmetries, const double mem_usage,
		  const double elapsed_time, const unsigned int times, const int procs)
{
  char kind_name[7], hostname[MAX_HOSTNAME_LENGTH];
  if(kind == APSP_NORMAL) strcpy(kind_name, "NORMAL");
  else                    strcpy(kind_name, "SAVING");
  gethostname(hostname, sizeof(hostname));
  time_t t = time(NULL);
  
  printf("------ Profile for APSP_RUN ------\n");
  printf("Date            = %s", ctime(&t));
  printf("Hostname        = %s\n", hostname);
  printf("Number of Times = %d\n", times);
  printf("Total Time      = %f sec.\n", elapsed_time);
  printf("Average Time    = %f sec.\n", elapsed_time/times);
  printf("Algorithm       = %s (%s)\n", kind_name, name);
  printf("Symmetries      = %d\n", symmetries);
  printf("Memory Usage    = %.3f MB\n", mem_usage);
  printf("Num of Procs    = %d\n", procs);
#ifdef _OPENMP
  printf("Num of Threads  = %d\n", omp_get_max_threads());
#else
  printf("Num of Threads  = %d\n", 1);
#endif
  printf("--------- End of Profile ---------\n");
}

bool apsp_check_loop(const int lines, int edge[lines][2])
{
  for(int i=0;i<lines;i++)
    if(edge[i][0] == edge[i][1]){
      printf("Loop is found in line %d\n", i+1);
      return false;
    }
  
  return true;
}

static bool has_multiple_edges(const int e00, const int e01, const int e10, const int e11)
{
  return ((e00 == e10 && e01 == e11) || (e00 == e11 && e01 == e10));
}

bool apsp_check_multiple_edges(const int lines, int edge[lines][2])
{
  for(int i=0;i<lines;i++)
    for(int j=i+1;j<lines;j++)
      if(has_multiple_edges(edge[i][0], edge[i][1], edge[j][0], edge[j][1])){
	printf("Multiple edeges are found in lines %d %d\n", i+1, j+1);
	return false;
      }
  
  return true;
}

int apsp_get_length(const int lines, const int edge[lines][2], const int height)
{
  int length = 0;
  for(int i=0;i<lines;i++)
    length = MAX(length, abs(edge[i][0]/height-edge[i][1]/height)+abs(edge[i][0]%height-edge[i][1]%height));
  
  return length;
}

bool apsp_check_general(char *fname)
{
  FILE *fp;
  if((fp = fopen(fname, "r")) == NULL)
    ERROR("File not found\n");

  int n1=1, n2=-1;
  fscanf(fp, "%d %d", &n1, &n2);
  fclose(fp);

  return (n2 != -1)? true : false;
}

int apsp_get_degree(const int nodes, const int lines, const int edge[lines][2])
{
  int node[nodes];
  for(int i=0;i<nodes;i++)
    node[i] = 0;

  for(int i=0;i<lines;i++){
    node[edge[i][0]]++;
    node[edge[i][1]]++;
  }

  int degree = node[0];
  for(int i=1;i<nodes;i++)
    degree = MAX(degree, node[i]);

  return degree;
}

int apsp_get_nodes(const int lines, const int (*edge)[2])
{
  int max = 0;
  for(int i=0;i<lines;i++){
    max = MAX(max, edge[i][0]);
    max = MAX(max, edge[i][1]);
  }

  return max + 1;
}

int apsp_get_lines(const char* fname)
{
  FILE *fp = NULL;
  if((fp = fopen(fname, "r")) == NULL)
    ERROR("File not found\n");

  int lines = 0, c;
  while((c = fgetc(fp)) != EOF)
    if(c == '\n')
      lines++;

  fclose(fp);

  return lines;
}

void apsp_read_edge_general(const char* fname, int (*edge)[2])
{
  FILE *fp;
  if((fp = fopen(fname, "r")) == NULL)
    ERROR("File not found\n");

  int n1, n2, i = 0;
  while(fscanf(fp, "%d %d", &n1, &n2) != EOF){
    edge[i][0] = n1;
    edge[i][1] = n2;
    i++;
  }

  fclose(fp);
}

void apsp_read_edge_grid(const char *fname, int *w, int *h, int (*edge)[2])
{
  FILE *fp;
  if((fp = fopen(fname, "r")) == NULL)
    ERROR("File not found\n");

  int n[4];
  *w = 0;
  *h = 0;
  while(fscanf(fp, "%d,%d %d,%d", &n[0], &n[1], &n[2], &n[3]) != EOF){
    *w = MAX(*w, n[0]);
    *h = MAX(*h, n[1]);
    *w = MAX(*w, n[2]);
    *h = MAX(*h, n[3]);
  }
  *w += 1;
  *h += 1;
  rewind(fp);

  int i = 0;
  while(fscanf(fp, "%d,%d %d,%d", &n[0], &n[1], &n[2], &n[3]) != EOF){
    edge[i][0] = n[0] * (*h) + n[1];
    edge[i][1] = n[2] * (*h) + n[3];
    i++;
  }

  fclose(fp);
}

// This function is inherited from "http://research.nii.ac.jp/graphgolf/py/create-random.py".
void apsp_set_lbounds_general(const int nodes, const int degree, int *low_diameter, double *low_ASPL)
{
  int diam = -1, n = 1, r = 1;
  double aspl = 0.0;

  while(1){
    int tmp = n + degree * pow(degree-1, r-1);
    if(tmp >= nodes)
      break;

    n = tmp;
    aspl += r * degree * pow(degree-1, r-1);
    diam = r++;
  }

  diam++;
  aspl += diam * (nodes - n);
  aspl /= (nodes - 1);

  *low_diameter = diam;
  *low_ASPL     = aspl;
}

static int dist(const int x1, const int y1, const int x2, const int y2)
{
  return(abs(x1 - x2) + abs(y1 - y2));
}

// This function is inherited from "http://research.nii.ac.jp/graphgolf/pl/lower-lattice.pl".
void apsp_set_lbounds_grid(const int m, const int n, const int degree, const int length, int *low_diameter, double *low_ASPL)
{
  int moore[m*n], hist[m*n], mh[m*n];
  int mn = m * n, current = degree, ii;
  double sum = 0;

  moore[0] = 1;
  moore[1] = degree + 1;
  for(ii=2;;ii++){
    current = current * (degree - 1);
    moore[ii] = moore[ii-1] + current;
    if(moore[ii] >= mn){
      moore[ii] = mn;
      break;
    }
  }

  int maxhop = MAX((m+n-2+(length-1))/length, ii);
  for(int i=ii+1;i<=maxhop;i++)
    moore[i] = mn;

  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      for(int k=0;k<=maxhop;k++)
        hist[k] = 0;

      for (int i2=0;i2<m;i2++)
        for(int j2=0;j2<n;j2++)
          hist[(dist(i,j,i2,j2)+length-1)/length]++;

      for(int k=1;k<=maxhop;k++)
        hist[k] += hist[k-1];

      for(int k=0;k<=maxhop;k++)
        mh[k] = MIN(hist[k], moore[k]);

      for(int k=1;k<=maxhop;k++)
        sum += (double)(mh[k] - mh[k-1]) * k;
    }
  }

  int dboth = 0;
  for(dboth=0;;dboth++)
    if(mh[dboth] == mn)
      break;

  *low_diameter = dboth;
  *low_ASPL     = sum/((double)mn*(mn-1));
}

void apsp_conv_adjacency2edge(const int nodes, const int degree, const int *num_degrees,
			      const int *adjacency, int (*edge)[2])
{
  char tmp[nodes][degree];
  for(int i=0;i<nodes;i++)
    for(int j=0;j<degree;j++)
      tmp[i][j] = NOT_VISITED;

  int j = 0;
  if(!num_degrees){
    for(int u=0;u<nodes;u++){
      for(int i=0;i<degree;i++){
        int v = *(adjacency + u * degree + i);
      	if(tmp[u][i] == NOT_VISITED){
          edge[j][0] = u;
          edge[j][1] = v;
	  tmp[u][i]  = VISITED;
          j++;
	  for(int k=0;k<degree;k++)
	    if(*(adjacency +v*degree+k) == u && tmp[v][k] == NOT_VISITED)
	      tmp[v][k] = VISITED;
        }
      }
    }
  }
  else{
    for(int u=0;u<nodes;u++){
      for(int i=0;i<num_degrees[u];i++){
	int v = *(adjacency + u * degree + i);
	if(tmp[u][i] == NOT_VISITED){
	  edge[j][0] = u;
	  edge[j][1] = v;
	  tmp[u][i]  = VISITED;
	  j++;
	  for(int k=0;k<num_degrees[v];k++)
            if(*(adjacency + v*degree+k) == u && tmp[v][k] == NOT_VISITED)
              tmp[v][k] = VISITED;
	}
      }
    }
  }
}

void apsp_conv_edge2adjacency(const int nodes, const int lines, const int edge[lines][2],
			      int *adjacency) // int adjacency[nodes][degree]
{
  int num_degrees[nodes];
  for(int i=0;i<nodes;i++)
    num_degrees[i] = 0;

  int degree = apsp_get_degree(nodes, lines, edge);
  for(int i=0;i<lines;i++){
    int n1 = edge[i][0];
    int n2 = edge[i][1];
    *(adjacency + n1 * degree + (num_degrees[n1]++)) = n2; //  adjacency[n1][num_degrees[n1]++] = n2;
    *(adjacency + n2 * degree + (num_degrees[n2]++)) = n1; //  adjacency[n2][num_degrees[n2]++] = n1;
  }
}

void apsp_conv_edge2adjacency_s(const int nodes, const int lines, const int edge[lines][2],
				const int symmetries, int *adjacency) // int adjacency[nodes/adjacency][degree]
{
  if(nodes % symmetries != 0)
    ERROR("nodes(%d) must be divisible by symmetries(%d)\n", nodes, symmetries);
  
  int based_nodes = nodes/symmetries;
  int num_degrees[based_nodes];
  for(int i=0;i<based_nodes;i++)
    num_degrees[i] = 0;

  int degree = apsp_get_degree(nodes, lines, edge);
  for(int i=0;i<lines;i++){
    int n1 = edge[i][0];
    int n2 = edge[i][1];
    if(n1 < based_nodes)
      *(adjacency + n1 * degree + (num_degrees[n1]++)) = n2; //  adjacency[n1][num_degrees[n1]++] = n2;
    if(n2 < based_nodes)
      *(adjacency + n2 * degree + (num_degrees[n2]++)) = n1; //  adjacency[n2][num_degrees[n2]++] = n1;
  }
}

void apsp_set_degrees(const int nodes, const int lines, int edge[lines][2],
		      int* num_degrees)
{
  for(int i=0;i<nodes;i++)
    num_degrees[i] = 0;

  for(int i=0;i<lines;i++){
    num_degrees[edge[i][0]]++;
    num_degrees[edge[i][1]]++;
  }
}

void apsp_mutate_adjacency_general(const int nodes, const int degree, const int *restrict num_degrees,
				   int adjacency[nodes][degree])
{
  int elements = 4; 
  int n[elements], d[elements];
  do{
    n[0] = get_random(nodes);
    n[1] = get_random(nodes);
  } while(n[0] == n[1]);
  
  do{
    d[0] = (!num_degrees)? get_random(degree) : get_random(num_degrees[n[0]]);
    n[2] = adjacency[n[0]][d[0]];
  } while(n[2] == n[1]);
  
  do{
    d[1] = (!num_degrees)? get_random(degree) : get_random(num_degrees[n[1]]);
    n[3] = adjacency[n[1]][d[1]];
  } while(n[3] == n[0]);

  if(n[0] == n[2]){ // loop
    for(int i=0;i<degree;i++)
      if(adjacency[n[2]][i] == n[0] && i != d[0])
	d[2] = i;
  }
  else{
    for(int i=0;i<degree;i++)
      if(adjacency[n[2]][i] == n[0])
	d[2] = i;
  }

  if(n[1] == n[3]){ // loop
    for(int i=0;i<degree;i++)
      if(adjacency[n[3]][i] == n[1] && i != d[1])
        d[3] = i;
  }
  else{
    for(int i=0;i<degree;i++)
      if(adjacency[n[3]][i] == n[1])
	d[3] = i;
  }

  if(!_n) _n = malloc(sizeof(int)*elements);
  if(!_d) _d = malloc(sizeof(int)*elements);
  memcpy(_n, n, sizeof(int)*elements);
  memcpy(_d, d, sizeof(int)*elements);
  _r = get_random(2);
  _nodes = nodes;
  _degree = degree;

  if(_r == 0){
    swap(&adjacency[n[0]][d[0]], &adjacency[n[1]][d[1]]);
    swap(&adjacency[n[2]][d[2]], &adjacency[n[3]][d[3]]);
  }
  else{
    swap(&adjacency[n[0]][d[0]], &adjacency[n[3]][d[3]]);
    swap(&adjacency[n[2]][d[2]], &adjacency[n[1]][d[1]]);
  }
}

void apsp_restore_adjacency(int adjacency[_nodes][_degree])
{
  if(_r == 0){
    swap(&adjacency[_n[0]][_d[0]], &adjacency[_n[1]][_d[1]]);
    swap(&adjacency[_n[2]][_d[2]], &adjacency[_n[3]][_d[3]]);
  }
  else{
    swap(&adjacency[_n[0]][_d[0]], &adjacency[_n[3]][_d[3]]);
    swap(&adjacency[_n[2]][_d[2]], &adjacency[_n[1]][_d[1]]);
  }
}

void apsp_generate_random_general(const int nodes, const int degree, int (*edge)[2])
{
  check_graph_parameters(nodes, degree);

  int half_degree = degree/2;
  for(int i=0;i<nodes-1;i++){
    for(int j=0;j<half_degree;j++){
      edge[i*half_degree+j][0] = i;
      edge[i*half_degree+j][1] = i+1;
    }
  }
  for(int j=0;j<half_degree;j++){
    int i = nodes - 1;
    edge[i*half_degree+j][0] = i;
    edge[i*half_degree+j][1] = 0;
  }

  if(degree%2 == 1){
    int half_node = nodes/2; // half_nodes must be a multiple of 2
    for(int i=0;i<half_node;i++){
      edge[half_degree*nodes+i][0] = i;
      edge[half_degree*nodes+i][1] = i+half_node;
    }
  }

  int lines = (nodes*degree)/2;
  int *adjacency = malloc(sizeof(int)*nodes*degree);
  apsp_conv_edge2adjacency(nodes, lines, edge, adjacency);

  for(int i=0;i<lines*GEN_GRAPH_ITERS;i++) // Give randomness
    apsp_mutate_adjacency_general(nodes, degree, NULL, (int (*)[degree])adjacency);

  apsp_conv_adjacency2edge(nodes, degree, NULL, adjacency, edge);

  free(adjacency);
}

void apsp_generate_random_general_s(const int nodes, const int degree, const int symmetries, int (*edge)[2])
{
  if(nodes % symmetries != 0)
    ERROR("nodes(%d) must be divisible by symmetries(%d)\n", nodes, symmetries);

  int based_nodes = nodes/symmetries;
  apsp_generate_random_general(based_nodes, degree, edge);
  int based_lines = (based_nodes * degree) / 2;

  for(int i=1;i<symmetries;i++){
    for(int j=0;j<based_lines;j++){
      for(int k=0;k<2;k++){
        int v = edge[j][k] + based_nodes * i;
        edge[i*based_lines+j][k] = (v < nodes)? v : v-nodes;
      }
    }
  }

  int *adjacency = malloc(sizeof(int) * based_nodes * degree);
  apsp_conv_edge2adjacency_s(nodes, (based_lines*symmetries), edge, symmetries, adjacency);
  // mutate_s x 100
  // adj 2 edge
  free(adjacency);
}
