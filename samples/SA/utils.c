#include "common.h"
static int _u[2], _v[2], _u_d[2], _v_d[2], _nodes, _degree, _symmetries, _kind, _rnd;

static void CHECK_PARAMETERS(const int nodes, const int degree)
{
  if(nodes % 2 == 1 && degree % 2 == 1)
    ERROR("Nodes(%d) or Degree(%d) must be a multiple of 2.\n", nodes, degree);
}

static int WIDTH(const int v, const int height)
{
  return  v/height;
}

static int HEIGHT(const int v, const int height)
{
  return v%height;
}

static int ROTATE(const int v, const int width, const int height, const int symmetries, const int degree)
{
  if(symmetries != 2 && symmetries != 4)
    ERROR("Invalid symmetries(%d)\n", symmetries);

  int w = WIDTH (v, height);
  int h = HEIGHT(v, height);
  if(symmetries == 2){
    if(degree != 180)
      ERROR("Invalid degree\n");

    return (width-w-1)*height + (height-h-1);
  }

  // symmetries == 4
  if(degree != 90 && degree != 180 && degree != 270)
    ERROR("Invalid degree\n");

  if(degree == 90)       return h*height + (height-w-1);
  else if(degree == 180) return (height-w-1)*height + (height-h-1);
  else                   return (height-h-1)*height + w; // degree == 270
}

static void SWAP(int *a, int *b)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

static bool IS_DIAMETER(const int u, const int v, const int nodes, const int symmetries)
{
  if(symmetries%2 != 0 || abs(u-v) != nodes/2)
    return false;
  else
    return true;
}

static void CHECK_SYMMETRIES(const int nodes, const int symmetries)
{
  if(nodes % symmetries != 0)
    ERROR("nodes(%d) must be divisible by symmetries(%d)\n", nodes, symmetries);
}

static void CHECK_SYMMETRIES_GRID(const int symmetries)
{
  if(symmetries != 1 && symmetries != 2 && symmetries != 4)
    ERROR("symmetries(%d) must be 1 or 2 or 4\n", symmetries);
}

static void CHECK_SYMMETRIES_WH(const int symmetries, const int width, const int height)
{
  if(symmetries == 2 && width%2 != 0)
    ERROR("width(%d) must be divisible by 2\n", width);
  else if(symmetries == 4 && (width%2 != 0 || height%2 != 0))
    ERROR("height(%d) must be divisible by 2\n", height);
  else if(symmetries == 4 && width != height)
    ERROR("Must be the same as width(%d) and height(%d)\n", width, height);
}

static int NORM(int x, const int nodes)
{
  while(x < 0 || x >= nodes)
    x = (x < 0)? x + nodes : x - nodes;
  
  return x;
}

// return adjacency[v][d];
static int GLOBAL_ADJ_GENERAL(const int nodes, const int degree, const int symmetries,
			      const int (*adjacency)[degree], const int v, const int d)
{
  int based_nodes = nodes/symmetries;
  int n = adjacency[v%based_nodes][d] + (v/based_nodes)*based_nodes;
  return NORM(n, nodes);
}

// Returns the local index of the vertices in the first quadrant
static int LOCAL_INDEX_GRID(const int x, const int height, const int symmetries)
{
  CHECK_SYMMETRIES_GRID(symmetries);
  
  if(symmetries == 2){
    return x;
  }
  else{ // symmetries == 4)
    int based_height = height/2;
    return WIDTH(x,height)*based_height + HEIGHT(x,height);
  }
}

// return adjacency[v][d];
static int GLOBAL_ADJ_GRID(const int width, const int height, const int degree, const int symmetries,
			   const int (*adjacency)[degree], const int v, const int d)
{
  CHECK_SYMMETRIES_GRID(symmetries);
  
  if(symmetries == 1){
    return adjacency[v][d];
  }
  else if(symmetries == 2){
    int based_width = width/2;
    if(WIDTH(v,height) < based_width)
      return adjacency[v][d];
    else{
      int w = adjacency[ROTATE(v, width, height, symmetries, 180)][d];
      return ROTATE(w, width, height, symmetries, 180);
    }
  }
  else{ // symmetries == 4
    int based_width  = width/2;
    int based_height = height/2;
    if(WIDTH(v,height) < based_width && HEIGHT(v,height) < based_height){
      return adjacency[LOCAL_INDEX_GRID(v,height,symmetries)][d];
    }
    else if(WIDTH(v,height) < based_width && HEIGHT(v,height) >= based_height){
      int x = ROTATE(v, width, height, symmetries, 270);
      int w = adjacency[LOCAL_INDEX_GRID(x,height,symmetries)][d];
      return ROTATE(w, width, height, symmetries, 90);
    }
    else if(WIDTH(v,height) >= based_width && HEIGHT(v,height) >= based_height){
      int x = ROTATE(v, width, height, symmetries, 180);
      int w = adjacency[LOCAL_INDEX_GRID(x,height,symmetries)][d];
      return ROTATE(w, width, height, symmetries, 180);
    }
    else{
      int x = ROTATE(v, width, height, symmetries, 90);
      int w = adjacency[LOCAL_INDEX_GRID(x,height,symmetries)][d];
      return ROTATE(w, width, height, symmetries, 270);
    }
  }
}

static int LOCAL_INDEX_general(const int v, const int position, const int nodes, const int symmetries)
{
  int based_nodes = nodes/symmetries;
  return NORM(v - (position/based_nodes)*based_nodes, nodes);
}

static int top_down_step(const int nodes, const int num_frontier, const int degree,
                         const int* restrict adjacency, int* restrict frontier,
                         int* restrict next, char* restrict bitmap)
{
  int count = 0;
  for(int i=0;i<num_frontier;i++){
    int v = frontier[i];
    for(int j=0;j<degree;j++){
      int n = *(adjacency + v * degree + j);  // int n = adjacency[v][j];
      if(bitmap[n] == NOT_VISITED){
        bitmap[n] = VISITED;
        next[count++] = n;
      }
    }
  }

  return count;
}

static bool simple_bfs(const int nodes, const int degree, int *adjacency)
{
  char *bitmap  = malloc(sizeof(char) * nodes);
  int *frontier = malloc(sizeof(int)  * nodes);
  int *next     = malloc(sizeof(int)  * nodes);
  int num_frontier = 1, root = 0;

  for(int i=0;i<nodes;i++)
    bitmap[i] = NOT_VISITED;

  frontier[0]  = root;
  bitmap[root] = VISITED;

  while(1){
    num_frontier = top_down_step(nodes, num_frontier, degree,
                                 adjacency, frontier, next, bitmap);
    if(num_frontier == 0) break;

    int *tmp = frontier;
    frontier = next;
    next     = tmp;
  }

  bool flag = false;  
  for(int i=0;i<nodes;i++)
    if(bitmap[i] == NOT_VISITED)
      flag = true;

  free(frontier);
  free(next);
  free(bitmap);
  
  return flag;
}

static int get_random(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
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

static bool has_multiple_edges(const int e00, const int e01, const int e10, const int e11)
{
  return ((e00 == e10 && e01 == e11) || (e00 == e11 && e01 == e10));
}

static int dist(const int x1, const int y1, const int x2, const int y2)
{
  return(abs(x1 - x2) + abs(y1 - y2));
}

static int get_lines(const int nodes, const int degree, const int *num_degrees)
{
  int lines = 0;
  if(!num_degrees){
    lines = (nodes * degree) / 2;
  }
  else{
    for(int i=0;i<nodes;i++)
      lines += num_degrees[i];
  }
  return lines;
}

static bool check_isolated_vertex(const int n, const int based_nodes, const int degree,
				  const int *num_degrees, const int (*adjacency)[degree])
{
  int x = n % based_nodes;
  int d = (!num_degrees)? degree : num_degrees[x];
  for(int i=1;i<d;i++)
    if(adjacency[x][0] != adjacency[x][i])
      return false;

  if(num_degrees)
    if(num_degrees[x] < num_degrees[adjacency[x][0]%based_nodes])
      return false; // It may not be an isolated vertex.
  
  return true;
}

void restore_adjacency(int *adjacency)
{
  int based_nodes = _nodes/_symmetries;
  adjacency[_u[0]%based_nodes * _degree + _u_d[0]] = LOCAL_INDEX_general(_v[0], _u[0], _nodes, _symmetries);
  adjacency[_v[0]%based_nodes * _degree + _v_d[0]] = LOCAL_INDEX_general(_u[0], _v[0], _nodes, _symmetries);
  if(_kind == MUTATE_1OPT) return;
  adjacency[_u[1]%based_nodes * _degree + _u_d[1]] = LOCAL_INDEX_general(_v[1], _u[1], _nodes, _symmetries);
  adjacency[_v[1]%based_nodes * _degree + _v_d[1]] = LOCAL_INDEX_general(_u[1], _v[1], _nodes, _symmetries);
}

static int get_degree_index_general(const int u, const int v, const int u_d, const int nodes,
				    const int symmetries, const int degree, const int (*adjacency)[degree])
{
  int based_nodes = nodes/symmetries;
  if(u == v){ // loop
    for(int i=0;i<degree;i++)
      if(GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, adjacency, v, i) == u && i != u_d)
	return i;
  }
  else if(symmetries%2 == 0 && abs(u-v) == nodes/2){
    return u_d;
  }
  else{
    for(int i=0;i<degree;i++)
      if(GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, adjacency, v, i) == u)
	return i;
  }

  ERROR("Something wrong ! [id=4]\n");
  return -1; // dummy
}

static int get_degree_index_grid(const int u, const int v, const int u_d, const int width, const int height,
				 const int symmetries, const int degree, const int (*adjacency)[degree])
{
  if(u == v){ // loop
    for(int i=0;i<degree;i++)
      if(GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, v, i) == u && i != u_d)
	return i;
  }
  else if(WIDTH(u,height) + WIDTH(v,height) == width-1 && HEIGHT(u,height) + HEIGHT(v,height) == height-1){
    return u_d;
  }
  else{
    for(int i=0;i<degree;i++){
      if(GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, v, i) == u)
	return i;
    }
  }
  
  ERROR("Something wrong ! [id=5]\n");
  return -1; // dummy
}


static void check_index(const int nodes, const int degree, const int symmetries,
			const int (*adjacency)[degree])
{
  int based_nodes = nodes/symmetries;
  for(int i=0;i<based_nodes;i++)
    for(int j=0;j<degree;j++)
      get_degree_index_general(i, adjacency[i][j], j, nodes, symmetries, degree, adjacency);
}

static void backup_restore_adjacency(const int u[2], const int u_d[2], const int v[2], const int v_d[2],
				     const int nodes, const int degree, const int symmetries, const int kind)
{
  int based_nodes = nodes/symmetries;
  memcpy(  _u,   u, sizeof(int)*2);
  memcpy(  _v,   v, sizeof(int)*2);
  memcpy(_u_d, u_d, sizeof(int)*2);
  memcpy(_v_d, v_d, sizeof(int)*2);
  _nodes      = nodes;
  _degree     = degree;
  _symmetries = symmetries;
  _kind       = kind;
}

bool mutate_adjacency_1opt_general_s(const int nodes, const int degree, const int *restrict num_degrees,
				     const int symmetries, int adjacency[nodes/symmetries][degree])
{
  int u[2], u_d[2], v[2], v_d[2]; // Declared with two elements since for backup_restore_adjacency()
  int based_nodes = nodes/symmetries;
  u[0]   = get_random(nodes);
  u_d[0] = (!num_degrees)? get_random(degree) : get_random(num_degrees[u[0]%based_nodes]);
  v[0]   = GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, adjacency, u[0], u_d[0]);
  if(symmetries%2 == 0 && abs(u[0]-v[0]) == nodes/2) return false;
  v_d[0] = get_degree_index_general(u[0], v[0], u_d[0], nodes, symmetries, degree, adjacency);
  backup_restore_adjacency(u, u_d, v, v_d, nodes, degree, symmetries, MUTATE_1OPT);

  // When it is an even number, there is one more pattern than when it is an odd number.
  // However, since the added pattern connects the vertices on the diagonal line,
  // As a result of some experiments, the pattern does not seem to be good, so comment it out.
  int rnd   = (symmetries%2 == 1)? get_random(symmetries-1) : get_random(symmetries);
  int new_v = (rnd != symmetries-1)? v[0] + based_nodes*(rnd+1) : u[0] + based_nodes*(symmetries/2);
  //  int rnd   = get_random(symmetries-1);
  //  int new_v = v[0] + based_nodes*(rnd+1);
  int tmp_v = adjacency[u[0]%based_nodes][u_d[0]];
  int new_u = v[0] - (new_v - u[0]);
  int tmp_u = adjacency[v[0]%based_nodes][v_d[0]];
  
  adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_INDEX_general(new_v, u[0], nodes, symmetries);
  adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_INDEX_general(new_u, v[0], nodes, symmetries);

  return true;
}

bool mutate_adjacency_2opt_general_s(const int nodes, const int degree, const int *restrict num_degrees,
				     const int symmetries, int adjacency[nodes/symmetries][degree])
{
  int u[2], v[2], u_d[2], v_d[2], based_nodes = nodes/symmetries;
  
  while(1){
    u[0] = get_random(nodes);
    u[1] = get_random(nodes);
    if(u[0] == u[1]) continue;
    
    u_d[0] = (!num_degrees)? get_random(degree) : get_random(num_degrees[u[0]%based_nodes]);
    v[0] = GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, adjacency, u[0], u_d[0]);
    if(v[0] == u[1]) continue;
    
    u_d[1] = (!num_degrees)? get_random(degree) : get_random(num_degrees[u[1]%based_nodes]);
    v[1] = GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, adjacency, u[1], u_d[1]);
    if(v[1] == u[0] || v[0] == v[1]) continue;
    break;
  }

  v_d[0] = get_degree_index_general(u[0], v[0], u_d[0], nodes, symmetries, degree, adjacency);
  v_d[1] = get_degree_index_general(u[1], v[1], u_d[1], nodes, symmetries, degree, adjacency);
  backup_restore_adjacency(u, u_d, v, v_d, nodes, degree, symmetries, MUTATE_2OPT);
  
  if(IS_DIAMETER(u[0], v[0], nodes, symmetries) && IS_DIAMETER(u[1], v[1], nodes, symmetries)){
    if((u[0] - u[1])%based_nodes == 0){
      return false;
    }
    else{
      adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_INDEX_general(u[1], u[0], nodes, symmetries);
      adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_INDEX_general(u[0], u[1], nodes, symmetries);
      return true;
    }
  }
  else if(IS_DIAMETER(u[0], v[0], nodes, symmetries) || IS_DIAMETER(u[1], v[1], nodes, symmetries)){
    if(IS_DIAMETER(u[1], v[1], nodes, symmetries)){
      SWAP(&u[0], &u[1]); SWAP(&u_d[0], &u_d[1]);
      SWAP(&v[0], &v[1]); SWAP(&v_d[0], &v_d[1]);
    }

    int opposite = nodes/2;
    int rnd = get_random(4);
    if(rnd == 0){
      adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_INDEX_general(v[1],          u[0], nodes, symmetries);
      adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_INDEX_general(u[1]+opposite, u[1], nodes, symmetries);
      adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_INDEX_general(v[1]+opposite, v[0], nodes, symmetries);
      adjacency[v[1]%based_nodes][v_d[1]] = LOCAL_INDEX_general(u[0],          v[1], nodes, symmetries);
    }
    else if(rnd == 1){
      adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_INDEX_general(v[1]+opposite, u[0], nodes, symmetries);
      adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_INDEX_general(u[1]+opposite, u[1], nodes, symmetries);
      adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_INDEX_general(v[1],          v[0], nodes, symmetries);
      adjacency[v[1]%based_nodes][v_d[1]] = LOCAL_INDEX_general(v[0],          v[1], nodes, symmetries);
    }
    else if(rnd == 2){
      adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_INDEX_general(u[1],          u[0], nodes, symmetries);
      adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_INDEX_general(u[0],          u[1], nodes, symmetries);
      adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_INDEX_general(u[1]+opposite, v[0], nodes, symmetries);
      adjacency[v[1]%based_nodes][v_d[1]] = LOCAL_INDEX_general(v[1]+opposite, v[1], nodes, symmetries);
    }
    else if(rnd == 3){
      adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_INDEX_general(u[1]+opposite, u[0], nodes, symmetries);
      adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_INDEX_general(v[0],          u[1], nodes, symmetries);
      adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_INDEX_general(u[1],          v[0], nodes, symmetries);
      adjacency[v[1]%based_nodes][v_d[1]] = LOCAL_INDEX_general(v[1]+opposite, v[1], nodes, symmetries);
    }
    return true;
  }
  else if((u[0]%based_nodes == u[1]%based_nodes && v[0]%based_nodes == v[1]%based_nodes) ||
	  (u[0]%based_nodes == v[1]%based_nodes && v[0]%based_nodes == u[1]%based_nodes)){
    // A graph with symmetry can be created when using mutate_adjacency_1opt_general_s().
    // But I want mutate_adjacency_1opt_general_s() not to be called in this function
    // because I want the number of calls to 1opt_s() and 2opt_s() to be about the same.
    return false;
  }
  _rnd = get_random(2);
  
  if(_rnd == 0){ // u[0]--v[1], v[0]--u[1]
    adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_INDEX_general(v[1], u[0], nodes, symmetries);
    adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_INDEX_general(v[0], u[1], nodes, symmetries);
    adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_INDEX_general(u[1], v[0], nodes, symmetries);
    adjacency[v[1]%based_nodes][v_d[1]] = LOCAL_INDEX_general(u[0], v[1], nodes, symmetries);
  }
  else{ // u[0]--u[1], v[0]--v[1]
    adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_INDEX_general(u[1], u[0], nodes, symmetries);
    adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_INDEX_general(u[0], u[1], nodes, symmetries);
    adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_INDEX_general(v[1], v[0], nodes, symmetries);
    adjacency[v[1]%based_nodes][v_d[1]] = LOCAL_INDEX_general(v[0], v[1], nodes, symmetries);
  }

  if(check_isolated_vertex(u[0], based_nodes, degree, num_degrees, adjacency) ||
     check_isolated_vertex(v[0], based_nodes, degree, num_degrees, adjacency) ||
     check_isolated_vertex(u[1], based_nodes, degree, num_degrees, adjacency) ||
     check_isolated_vertex(v[1], based_nodes, degree, num_degrees, adjacency)){
    restore_adjacency((int *)adjacency);
    return false;
  }
  else{
    return true;
  }
}

void mutate_adjacency_general_s(const int nodes, const int degree, const int *restrict num_degrees,
				const int symmetries, int adjacency[nodes/symmetries][degree])
{
  CHECK_SYMMETRIES(nodes, symmetries);
  if(symmetries == 1){
    while(1){
      if(mutate_adjacency_2opt_general_s(nodes, degree, num_degrees, symmetries, adjacency))
	break;
    }
    return;
  }
  
  while(1){
    if(get_random(2) == 0){
      if(mutate_adjacency_2opt_general_s(nodes, degree, num_degrees, symmetries, adjacency))
	break;
    }
    else{
      if(mutate_adjacency_1opt_general_s(nodes, degree, num_degrees, symmetries, adjacency))
	break;
    }
  }
}

void mutate_adjacency_general(const int nodes, const int degree, const int *restrict num_degrees,
			      int adjacency[nodes][degree])
{
  mutate_adjacency_general_s(nodes, degree, num_degrees, 1, adjacency);
}

void mutate_adjacency_grid(const int width, const int height, const int degree,
			   const int *restrict num_degrees, const int length, int (*adjacency)[degree])
{
  int nodes = width * height;
  while(1){
    mutate_adjacency_general(nodes, degree, num_degrees, adjacency);
    if(_rnd == 0 && check_length(_u[0], _v[1], height, length) && check_length(_u[1], _v[0], height, length)){
      break;
    }
    else if(_rnd == 1 && check_length(_u[0], _u[1], height, length) && check_length(_v[0], _v[1], height, length)){
      break;
    }
    else{
      restore_adjacency((int *)adjacency);
    }
  }
}

static bool mutate_adjacency_1opt_grid_s(const int width, const int height, const int degree,
					 const int *restrict num_degrees, const int length, const int symmetries,
					 int (*adjacency)[degree])
{
  int nodes       = width * height;
  int based_nodes = nodes/symmetries;
  int u           = get_random(nodes);
  int u_d         = (!num_degrees)? get_random(degree) : get_random(num_degrees[u%based_nodes]);
  int v           = GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, u, u_d);
  if(WIDTH(u,height) + WIDTH(v,height) == width-1 && HEIGHT(u,height) + HEIGHT(v,height) == height-1)
    return false;
  int v_d         = get_degree_index_grid(u, v, u_d, width, height, symmetries, degree, adjacency);

  int rnd   = get_random(symmetries-1);
  int new_v = ROTATE(v, width, height, symmetries, (rnd+1)*(360/symmetries));
  if(!check_length(u, new_v, height, length)) return false;
  int new_u = ROTATE(u, width, height, symmetries, (symmetries-rnd-1)*(360/symmetries));
  int based_width  = width/2;
  int based_height = (symmetries == 2)? height : height/2;

  if(symmetries == 2){
    // Update v
    if(WIDTH(new_v,height) < based_width){
      adjacency[LOCAL_INDEX_GRID(new_v,height,symmetries)][v_d] = u;
    }
    else{
      int x = ROTATE(new_v, width, height, symmetries, 180);
      adjacency[LOCAL_INDEX_GRID(x,height,symmetries)][v_d] = ROTATE(u, width, height, symmetries, 180);
    }

    // Update u
    if(WIDTH(new_u,height) < based_width){
      adjacency[LOCAL_INDEX_GRID(new_u,height,symmetries)][u_d] = v;
    }
    else{
      int x = ROTATE(new_u, width, height, symmetries, 180);
      adjacency[LOCAL_INDEX_GRID(x,height,symmetries)][u_d] = ROTATE(v, width, height, symmetries, 180);
    }
  }
  else{ // symmetries == 4
    // Update v
    if(WIDTH(new_v,height) < based_width && HEIGHT(new_v,height) < based_height){
      adjacency[LOCAL_INDEX_GRID(new_v,height,symmetries)][v_d] = u;
    }
    else if(WIDTH(new_v,height) < based_width && HEIGHT(new_v,height) >= based_height){
      int x = ROTATE(new_v, width, height, symmetries, 270);
      adjacency[LOCAL_INDEX_GRID(x,height,symmetries)][v_d] = ROTATE(u, width, height, symmetries, 270);
    }
    else if(WIDTH(new_v,height) >= based_width && HEIGHT(new_v,height) >= based_height){
      int x = ROTATE(new_v, width, height, symmetries, 180);
      adjacency[LOCAL_INDEX_GRID(x,height,symmetries)][v_d] = ROTATE(u, width, height, symmetries, 180);
    }
    else{
      int x = ROTATE(new_v, width, height, symmetries, 90);
      adjacency[LOCAL_INDEX_GRID(x,height,symmetries)][v_d] = ROTATE(u, width, height, symmetries,  90);
    }

    // Update u
    if(WIDTH(new_u,height) < based_width && HEIGHT(new_u,height) < based_height){
      adjacency[LOCAL_INDEX_GRID(new_u,height,symmetries)][u_d] = v;
    }
    else if(WIDTH(new_u,height) < based_width && HEIGHT(new_u,height) >= based_height){
      int x = ROTATE(new_u, width, height, symmetries, 270);
      adjacency[LOCAL_INDEX_GRID(x,height,symmetries)][u_d] = ROTATE(v, width, height, symmetries, 270);
    }
    else if(WIDTH(new_u,height) >= based_width && HEIGHT(new_u,height) >= based_height){
      int x = ROTATE(new_u, width, height, symmetries, 180);
      adjacency[LOCAL_INDEX_GRID(x,height,symmetries)][u_d] = ROTATE(v, width, height, symmetries, 180);
    }
    else{
      int x = ROTATE(new_u, width, height, symmetries, 90);
      adjacency[LOCAL_INDEX_GRID(x,height,symmetries)][u_d] = ROTATE(v, width, height, symmetries,  90);
    }
  }
  return true;
}

void mutate_adjacency_grid_s(const int width, const int height, const int degree,
			     const int *restrict num_degrees, const int length, const int symmetries,
			     int (*adjacency)[degree])
{
  CHECK_PARAMETERS(width*height, degree);
  CHECK_SYMMETRIES_GRID(symmetries);  
  CHECK_SYMMETRIES_WH(symmetries, width, height);
  
  if(symmetries == 1){
    mutate_adjacency_grid(width, height, degree, num_degrees, length, adjacency);
    return;
  }

  while(1){
    if(mutate_adjacency_1opt_grid_s(width, height, degree, num_degrees, length, symmetries, adjacency))
      break;
  }
}
