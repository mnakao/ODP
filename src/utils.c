#include "common.h"
static time_t _start_t;
static double _elapsed_time;

static void check_graph_parameters(const int nodes, const int degree)
{
  if(nodes % 2 == 1 && degree % 2 == 1)
    ERROR("Nodes(%d) or Degree(%d) must be a multiple of 2.\n", nodes, degree);
}

static int get_random(const int max)
{
  return (int)(random()*((double)max)/(1.0+RAND_MAX));
}

static void swap(int *a, int *b)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

static bool check_duplicated_vertex(const int e00, const int e01, const int e10, const int e11)
{
  return (e00 == e10 || e01 == e11 || e00 == e11 || e01 == e10);
}

static void simple_2opt_general(const int lines, int edge[lines][2])
{
  int e0, e1;
  do{
    e0 = get_random(lines);
    e1 = get_random(lines);
  } while(check_duplicated_vertex(edge[e0][0], edge[e0][1], edge[e1][0], edge[e1][1]));

  if(get_random(2) == 0)
    swap(&edge[e0][1], &edge[e1][1]);
  else
    swap(&edge[e0][1], &edge[e1][0]);
}

void apsp_random_general(const int nodes, const int degree, const unsigned int seed, int *edge)
{
  check_graph_parameters(nodes, degree);
  srand(seed);

  int half_degree = degree/2;
  for(int i=0;i<nodes-1;i++){
    for(int j=0;j<half_degree;j++){
      edge[(i*half_degree+j)*2  ] = i;
      edge[(i*half_degree+j)*2+1] = i+1;
    }
  }
  for(int j=0;j<half_degree;j++){
    int i = nodes - 1;
    edge[(i*half_degree+j)*2  ] = i;
    edge[(i*half_degree+j)*2+1] = 0;
  }

  if(degree%2 == 1){
    int half_node = nodes/2; // half_nodes must be a multiple of 2
    for(int i=0;i<half_node;i++){
      edge[(half_degree*nodes+i)*2  ] = i;
      edge[(half_degree*nodes+i)*2+1] = i+half_node;
    }
  }

  int lines = (nodes*degree)/2;
  for(int i=0;i<lines*GEN_GRAPH_ITERS;i++) // Give randomness
    simple_2opt_general(lines, (int (*)[2])edge);
}

void apsp_output_edge_general(char *fname, const int lines, const int edge[lines][2])
{
  FILE *fp = NULL;
  
  if((fp = fopen(fname, "w")) == NULL)
    ERROR("Cannot open %s\n", fname);
  
  for(int i=0;i<lines;i++)
    fprintf(fp, "%d %d\n", edge[i][0], edge[i][1]);
  
  fclose(fp);
}

void apsp_output_edge_grid(char *fname, const int lines, const int height, const int edge[lines][2])
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

void apsp_random_general_s(const int nodes, const int degree, const int groups, const unsigned int seed, int *edge)
{
  check_graph_parameters(nodes, degree);
  srand(seed);
  ERROR("Not implemented yet\n");
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
    } while(check_duplicated_vertex(edge[e0][0], edge[e0][1], edge[e1][0], edge[e1][1]));

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
void apsp_random_grid(const int width, const int height, const int degree,
		      const int length, const unsigned int seed, int *edge)
{
  int nodes = width * height;
  check_graph_parameters(nodes, degree);
  srand(seed);

  int i = 0;
  for(int x=0;x<width/2;x++){
    for(int y=0;y<height;y++){
      for(int k=0;k<degree;k++){
        edge[i*2]   = y + 2 * x * height;
        edge[i*2+1] = edge[i*2] + height;
        i++;
      }
    }
  }

  if(width%2 == 1){
    for(int y=0;y<height/2;y++){
      for(int k=0;k<degree;k++){
        edge[i*2]   = (width - 1) * height + 2 * y;
        edge[i*2+1] = edge[i*2] + 1;
        i++;
      }
    }

    /* add self-loop */
    if(height%2 == 1){
      for(int k=0;k<degree/2;k++){
        edge[i*2] = edge[i*2+1] = nodes - 1;
        i++;
      }
    }
  }

  int lines = (nodes*degree)/2;
  for(int i=0;i<lines*GEN_GRAPH_ITERS;i++)  // Give randomness
    simple_2opt_grid(height, length, lines, (int (*)[2])edge);
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

void matmul(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
	    const int *num_degrees, const int *restrict adjacency, const bool enable_avx2, const int elements)
{
  if(enable_avx2){
    if(!num_degrees){
#pragma omp parallel for
      for(int i=0;i<nodes;i++){
	__m256i *b = (__m256i *)(B + i*elements);
	for(int j=0;j<degree;j++){
	  int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
	  __m256i *a = (__m256i *)(A + n*elements);
	  for(int k=0;k<elements/4;k++){
	    __m256i aa = _mm256_load_si256(a+k);
	    __m256i bb = _mm256_load_si256(b+k);
	    _mm256_store_si256(b+k, _mm256_or_si256(aa, bb));
	  }
	}
      }
    }
    else{
#pragma omp parallel for
      for(int i=0;i<nodes;i++){
	__m256i *b = (__m256i *)(B + i*elements);
	for(int j=0;j<num_degrees[i];j++){
	  int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
	  __m256i *a = (__m256i *)(A + n*elements);
	  for(int k=0;k<elements/4;k++){
	    __m256i aa = _mm256_load_si256(a+k);
	    __m256i bb = _mm256_load_si256(b+k);
	    _mm256_store_si256(b+k, _mm256_or_si256(aa, bb));
	  }
	}
      }
    }
  }
  else{
    if(!num_degrees){
#pragma omp parallel for
      for(int i=0;i<nodes;i++){
	for(int j=0;j<degree;j++){
	  int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
	  for(int k=0;k<elements;k++)
	    B[i*elements+k] |= A[n*elements+k];
	}
      }
    }
    else{
#pragma omp parallel for
      for(int i=0;i<nodes;i++){
	for(int j=0;j<num_degrees[i];j++){
	  int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
	  for(int k=0;k<elements;k++)
	    B[i*elements+k] |= A[n*elements+k];
	}
      }
    }
  }
}

void matmul_CHUNK(const uint64_t *restrict A, uint64_t *restrict B, const int nodes, const int degree,
		  const int *num_degrees, const int *restrict adjacency, const bool enable_avx2)
{
  if(enable_avx2){
    if(!num_degrees){
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
    else{
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
  }
  else{
    if(!num_degrees){
#pragma omp parallel for
      for(int i=0;i<nodes;i++){
	for(int j=0;j<degree;j++){
	  int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
	  for(int k=0;k<CPU_CHUNK;k++)
	    B[i*CPU_CHUNK+k] |= A[n*CPU_CHUNK+k];
	}
      }
    }
    else{
#pragma omp parallel for
      for(int i=0;i<nodes;i++){
	for(int j=0;j<num_degrees[i];j++){
	  int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
	  for(int k=0;k<CPU_CHUNK;k++)
	    B[i*CPU_CHUNK+k] |= A[n*CPU_CHUNK+k];
	}
      }
    }
  }
}

double apsp_get_mem_usage(const int kind, const int nodes, const int degree, const int groups,
			  const int *num_degree, const int procs, const int chunk)
{
  double edge_size   = (double)nodes * degree * sizeof(int);
  double adj_size    = edge_size;
  double degree_size = (num_degree)? (double)nodes * sizeof(int) : 0;
  double normal_mem  = (double)nodes * ((double)nodes / (4 * groups * procs));
  double saving_mem  = 16 * (double)nodes * chunk;
  double in_apsp     = (kind == APSP_NORMAL)? normal_mem : saving_mem;
  
  return (edge_size + adj_size + degree_size + in_apsp)/1024/1024;
}

int apsp_get_kind(const int nodes, const int degree, const int* num_degrees, const int groups,
		  const int procs, const int chunk)
{
  char *val = getenv("APSP");
  int kind;
  if(val == NULL){
    double normal_mem_usage = apsp_get_mem_usage(APSP_NORMAL, nodes, degree, groups, num_degrees, procs, chunk);
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

static bool check_profile()
{
  char *val = getenv("APSP_PROFILE");
  if(val){
    if(atoi(val) == 1)
      return true;
  }
  return false;
}

static double get_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1.0e-6 * t.tv_usec;
}


void apsp_start_profile()
{
  if(!check_profile()) return;
  
  _elapsed_time = get_time();
  _start_t      = time(NULL);
}

void apsp_end_profile(const char* name, const int kind, const int groups, const double mem_usage, const int procs)
{
  if(!check_profile()) return;
  
  time_t end_t = time(NULL);
  char kind_name[7], hostname[MAX_HOSTNAME_LENGTH];
  if(kind == APSP_NORMAL) strcpy(kind_name, "NORMAL");
  else                    strcpy(kind_name, "SAVING");
  gethostname(hostname, sizeof(hostname));
  
  printf("------ Start of Profile ------\n");
  printf("Hostname       = %s\n", hostname);
  printf("Start Date     = %s", ctime(&_start_t));
  printf("End Date       = %s", ctime(&end_t));
  printf("Elapsed Time   = %f sec.\n", get_time() - _elapsed_time);
  printf("Algorithm      = %s (%s)\n", kind_name, name);
  printf("Symmetry       = %d\n", groups);
  printf("Memory Usage   = %.3f MB\n", mem_usage);
  if(procs != 1) printf("Procs          = %d\n", procs);
#ifdef _OPENMP
  printf("Num of Threads = %d\n", omp_get_max_threads());
#endif
  printf("------  End of Profile  ------\n");
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

static bool has_duplicated_edge(const int e00, const int e01, const int e10, const int e11)
{
  return ((e00 == e10 && e01 == e11) || (e00 == e11 && e01 == e10));
}

bool apsp_check_duplicated_edge(const int lines, int edge[lines][2])
{
  for(int i=0;i<lines;i++)
    for(int j=i+1;j<lines;j++)
      if(has_duplicated_edge(edge[i][0], edge[i][1], edge[j][0], edge[j][1])){
	printf("Duplicate edeges are found in lines %d %d\n", i+1, j+1);
	return false;
      }
  
  return true;
}

int apsp_get_length(const int lines, int edge[lines][2], const int height)
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

int apsp_get_degree(const int nodes, const int lines, int edge[lines][2])
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

int apsp_get_nodes(const int lines, int (*edge)[2])
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

void apsp_set_edge_general(const char* fname, int (*edge)[2])
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

void apsp_set_edge_grid(const char *fname, int *w, int *h, int (*edge)[2])
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

void apsp_set_adjacency(const int nodes, const int degree, const int lines,
			const int edge[lines][2], int *adjacency)
{
  int num_degrees[nodes];
  for(int i=0;i<nodes;i++)
    num_degrees[i] = 0;
  
  for(int i=0;i<lines;i++){
    int n1 = edge[i][0];
    int n2 = edge[i][1];
    *(adjacency + n1 * degree + (num_degrees[n1]++)) = n2; //  adjacency[n1][num_degrees[n1]++] = n2;
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
