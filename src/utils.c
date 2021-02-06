#include "common.h"
int _u[2], _v[2], _u_d[2], _v_d[2];

static int WIDTH(const int v, const int height)
{
  return  v/height;
}

static int HEIGHT(const int v, const int height)
{
  return v%height;
}

void ODP_Print_adjacency(const int nodes, const int degree, const int num_degrees[nodes], const int adjacency[nodes][degree])
{
  for(int i=0;i<nodes;i++){
    printf("%d ", i);
    int d = (!num_degrees)? degree : num_degrees[i];
    for(int j=0;j<d;j++){
      if(adjacency[i][j] != NOT_DEFINED)
	printf("%3d", adjacency[i][j]);
    }
    printf("\n");
  }
}

void ODP_Print_edge_general(const int lines, const int edge[lines][2])
{
  for(int i=0;i<lines;i++)
    printf("%d %d\n", edge[i][0], edge[i][1]);
}

void ODP_Print_edge_grid(const int lines, const int height, const int edge[lines][2])
{
  for(int i=0;i<lines;i++)
    printf("%d,%d %d,%d\n",
           WIDTH(edge[i][0],height), HEIGHT(edge[i][0],height),
	   WIDTH(edge[i][1],height), HEIGHT(edge[i][1],height));
}

static int DISTANCE(const int u, const int v, const int height)
{
  int u_w = WIDTH(u,height);
  int u_h = HEIGHT(u,height);
  int v_w = WIDTH(v,height);
  int v_h = HEIGHT(v,height);
  return abs(u_w - v_w) + abs(u_h - v_h);
}

static bool CHECK_LENGTH(const int u, const int v, const int height, const int length)
{

  return (DISTANCE(u, v, height) <= length);
}

int ODP_ROTATE(const int v, const int width, const int height, const int symmetries, const int degree)
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

static bool IS_DIAMETER_GENERAL(const int u, const int v, const int nodes, const int symmetries)
{
  if(symmetries%2 != 0 || abs(u-v) != nodes/2)
    return false;
  else
    return true;
}

static bool IS_DIAMETER_GRID(const int u, const int v, const int width, const int height, const int symmetries)
{
  if(symmetries == 1)
    return false;
  else
    return (u == ODP_ROTATE(v, width, height, symmetries, 180));
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

double ODP_Get_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1.0e-6 * t.tv_usec;
}

bool ODP_Check_profile()
{
  char *val = getenv("ODP_PROFILE");
  if(val){
    if(atoi(val) == 1)
      return true;
  }
  return false;
}

void ODP_Conv_adjacency2edge_grid_s(const int width, const int height, const int degree, const int *num_degrees,
                                    const int (*adjacency)[degree], const int symmetries, int (*edge)[2])
{
  CHECK_SYMMETRIES_GRID(symmetries);
  CHECK_SYMMETRIES_WH(symmetries, width, height);

  int nodes = width * height;
  if(symmetries == 1){
    ODP_Conv_adjacency2edge(nodes, degree, num_degrees, (int *)adjacency, edge);
    return;
  }

  int (*tmp)[degree] = malloc(sizeof(int) * nodes * degree);
  int based_nodes = nodes/symmetries;
  if(symmetries == 2){
    for(int i=0;i<based_nodes;i++){
      int d = (!num_degrees)? degree : num_degrees[i];
      for(int k=0;k<d;k++){
        int v = adjacency[i][k];
        tmp[i][k] = v;
        tmp[ODP_ROTATE(i, width, height, symmetries, 180)][k] = ODP_ROTATE(v, width, height, symmetries, 180);
      }
    }
  }
  else{ // symmetries == 4
    int based_height = height/2;
    for(int i=0;i<based_nodes;i++){
      int d = (!num_degrees)? degree : num_degrees[i];
      int j = (i/based_height)*height + (i%based_height);
      for(int k=0;k<d;k++){
        int v = adjacency[i][k];
        tmp[j][k] = v;
        tmp[ODP_ROTATE(j, width, height, symmetries,  90)][k] = ODP_ROTATE(v, width, height, symmetries,  90);
        tmp[ODP_ROTATE(j, width, height, symmetries, 180)][k] = ODP_ROTATE(v, width, height, symmetries, 180);
        tmp[ODP_ROTATE(j, width, height, symmetries, 270)][k] = ODP_ROTATE(v, width, height, symmetries, 270);
      }
    }
  }
  ODP_Conv_adjacency2edge(nodes, degree, num_degrees, (int *)tmp, edge);

  free(tmp);
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
int ODP_LOCAL_INDEX_GRID(const int x, const int width, const int height, const int symmetries)
{
  CHECK_SYMMETRIES_GRID(symmetries);

  int based_width = width/2;
  if(symmetries == 1){
    return x;
  }
  else if(symmetries == 2){
    if(WIDTH(x,height) >= based_width)
      ERROR("Something Wrong ! [id=6] %d %d %d\n", x, height, based_width);
    
    return x;
  }
  else{ // symmetries == 4
    int based_height = height/2;
    if(WIDTH(x,height) >= based_width || HEIGHT(x,height) >= based_height)
       ERROR("Something Wrong ! [id=7]\n");
    
    return WIDTH(x,height)*based_height + HEIGHT(x,height);
  }
}

// return adjacency[v][d];
int ODP_GLOBAL_ADJ_GRID(const int width, const int height, const int degree, const int symmetries,
			const int (*adjacency)[degree], const int v, const int d)
{
  if(symmetries == 1){
    return adjacency[v][d];
  }
  else if(symmetries == 2){
    int based_width = width/2;
    if(WIDTH(v,height) < based_width)
      return adjacency[v][d];
    else{
      int y = adjacency[ODP_ROTATE(v, width, height, symmetries, 180)][d];
      return ODP_ROTATE(y, width, height, symmetries, 180);
    }
  }
  else{ // symmetries == 4
    int based_width  = width/2;
    int based_height = height/2;
    if(WIDTH(v,height) < based_width && HEIGHT(v,height) < based_height){
      return adjacency[ODP_LOCAL_INDEX_GRID(v,width,height,symmetries)][d];
    }
    else if(WIDTH(v,height) < based_width && HEIGHT(v,height) >= based_height){
      int x = ODP_ROTATE(v, width, height, symmetries, 270);
      int y = adjacency[ODP_LOCAL_INDEX_GRID(x,width,height,symmetries)][d];
      return ODP_ROTATE(y, width, height, symmetries, 90);
    }
    else if(WIDTH(v,height) >= based_width && HEIGHT(v,height) >= based_height){
      int x = ODP_ROTATE(v, width, height, symmetries, 180);
      int y = adjacency[ODP_LOCAL_INDEX_GRID(x,width,height,symmetries)][d];
      return ODP_ROTATE(y, width, height, symmetries, 180);
    }
    else{
      int x = ODP_ROTATE(v, width, height, symmetries, 90);
      int y = adjacency[ODP_LOCAL_INDEX_GRID(x,width,height,symmetries)][d];
      return ODP_ROTATE(y, width, height, symmetries, 270);
    }
  }
}

static int LOCAL_INDEX_GENERAL(const int v, const int position, const int nodes, const int symmetries)
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

static void CHECK_PARAMETERS(const int nodes, const int degree)
{
  if(nodes % 2 == 1 && degree % 2 == 1)
    ERROR("Nodes(%d) or Degree(%d) must be a multiple of 2.\n", nodes, degree);
}

static int get_random(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
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

static int global2local_vertex_grid(const int v, const int width, const int height, const int symmetries)
{
  CHECK_SYMMETRIES_GRID(symmetries);
  
  if(symmetries == 1){
    return v;
  }
  else if(symmetries == 2){
    int based_width = width/2;
    if(WIDTH(v,height) < based_width)
      return v;
    else
      return ODP_ROTATE(v,width,height,symmetries,180);
  }
  else{ //  symmetries == 4
    int based_width  = width/2;
    int based_height = height/2;
    if(WIDTH(v,height) < based_width && HEIGHT(v,height) < based_height){
      return ODP_LOCAL_INDEX_GRID(v,width,height,symmetries);
    }
    else if(WIDTH(v,height) < based_width && HEIGHT(v,height) >= based_height){
      return ODP_LOCAL_INDEX_GRID(ODP_ROTATE(v,width,height,symmetries,270),width,height,symmetries);
    }
    else if(WIDTH(v,height) >= based_width && HEIGHT(v,height) >= based_height){
      return ODP_LOCAL_INDEX_GRID(ODP_ROTATE(v,width,height,symmetries,180),width,height,symmetries);
    }
    else{
      return ODP_LOCAL_INDEX_GRID(ODP_ROTATE(v,width,height,symmetries,90),width,height,symmetries);
    }
  }
}

static bool check_isolated_vertex_grid(const int v, const int width, const int height, const int symmetries,
				       const int degree, const int *num_degrees, const int (*adjacency)[degree])
{
  int based_width = width/2;
  int based_height = (symmetries == 2)? height : height/2;
  int x = global2local_vertex_grid(v, width, height, symmetries);
  int d = (!num_degrees)? degree : num_degrees[x]; // x is a local vertex
  for(int i=1;i<d;i++)
    if(adjacency[x][0] != adjacency[x][i])
      return false;

  if(num_degrees){
    if(num_degrees[x] < num_degrees[global2local_vertex_grid(adjacency[x][0], width, height, symmetries)])
      return false;
  }

  return true;
}

static bool check_isolated_vertex_general(const int n, const int based_nodes, const int degree,
					  const int *num_degrees, const int (*adjacency)[degree])
{
  int x = n % based_nodes;
  int d = (!num_degrees)? degree : num_degrees[x];
  for(int i=1;i<d;i++)
    if(adjacency[x][0] != adjacency[x][i])
      return false;

  if(num_degrees)
    if(num_degrees[x] < num_degrees[adjacency[x][0]%based_nodes])
      return false;
  
  return true;
}

// adjacency[v][v_d] = u;
static void set_adjacency(const int v, const int v_d, const int u, const int width, const int height,
                          const int degree, const int symmetries, int (*adjacency)[degree])
{
  int based_width  = width/2;
  int based_height = (symmetries == 2)? height : height/2;
  if(symmetries == 1){
    adjacency[v][v_d] = u;
  }
  else if(symmetries == 2){
    if(WIDTH(v,height) < based_width){
      adjacency[v][v_d] = u;
    }
    else{
      int x = ODP_ROTATE(v, width, height, symmetries, 180);
      adjacency[x][v_d] = ODP_ROTATE(u, width, height, symmetries, 180);
    }
  }
  else if(symmetries == 4){
    if(WIDTH(v,height) < based_width && HEIGHT(v,height) < based_height){
      adjacency[ODP_LOCAL_INDEX_GRID(v,width,height,symmetries)][v_d] = u;
    }
    else if(WIDTH(v,height) < based_width && HEIGHT(v,height) >= based_height){
      int x = ODP_ROTATE(v, width, height, symmetries, 270);
      adjacency[ODP_LOCAL_INDEX_GRID(x,width,height,symmetries)][v_d] = ODP_ROTATE(u, width, height, symmetries, 270);
    }
    else if(WIDTH(v,height) >= based_width && HEIGHT(v,height) >= based_height){
      int x = ODP_ROTATE(v, width, height, symmetries, 180);
      adjacency[ODP_LOCAL_INDEX_GRID(x,width,height,symmetries)][v_d] = ODP_ROTATE(u, width, height, symmetries, 180);
    }
    else{
      int x = ODP_ROTATE(v, width, height, symmetries, 90);
      adjacency[ODP_LOCAL_INDEX_GRID(x,width,height,symmetries)][v_d] = ODP_ROTATE(u, width, height, symmetries,  90);
    }
  }
}

void ODP_Restore_adjacency_grid_s(const int width, const int height, const int degree,
				  const int symmetries, int (*adjacency)[degree])
{
  for(int i=0;i<2;i++){
    set_adjacency(_u[i], _u_d[i], _v[i], width, height, degree, symmetries, adjacency);
    set_adjacency(_v[i], _v_d[i], _u[i], width, height, degree, symmetries, adjacency);
  }
}

void ODP_Restore_adjacency_grid(const int width, const int height, const int degree,
				int (*adjacency)[degree])
{
  ODP_Restore_adjacency_grid_s(width, height, degree, 1, adjacency);
}

static bool check_multiple_edges_general_s(const int u, const int u_d, const int v, const int nodes, const int degree,
					   const int *num_degrees, const int symmetries, const int (*adjacency)[degree])
{
  int based_nodes = nodes/symmetries;
  int d = (!num_degrees)? degree : num_degrees[u%based_nodes];
  for(int i=0;i<d;i++)
    if(i!=u_d && adjacency[u][i] == v)
      return false;

  return true;
}

static bool check_multiple_edges_grid_s(const int u, const int u_d, const int v, const int width, const int height,
					const int degree, const int *num_degrees, const int symmetries, const int (*adjacency)[degree])
{
  int d = (!num_degrees)? degree : num_degrees[global2local_vertex_grid(u, width, height, symmetries)];
  for(int i=0;i<d;i++)
    if(i!=u_d && ODP_GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, u, i) == v)
      return false;

  return true;
}

static bool check_isolated_vertices_grid(const int u[2], const int v[2], const int width, const int height, const int symmetries,
                                         const int degree, const int *num_degrees, int (*adjacency)[degree])
{
  if(check_isolated_vertex_grid(u[0], width, height, symmetries, degree, num_degrees, adjacency) ||
     check_isolated_vertex_grid(v[0], width, height, symmetries, degree, num_degrees, adjacency) ||
     check_isolated_vertex_grid(u[1], width, height, symmetries, degree, num_degrees, adjacency) ||
     check_isolated_vertex_grid(v[1], width, height, symmetries, degree, num_degrees, adjacency)){
    ODP_Restore_adjacency_grid_s(width, height, degree, symmetries, adjacency);
    return false;
  }
  else
    return true;
}

void ODP_Restore_adjacency_general_s(const int nodes, const int degree, const int symmetries, 
				     int (*adjacency)[degree])
{
  int based_nodes = nodes/symmetries;
  for(int i=0;i<2;i++){
    adjacency[_u[i]%based_nodes][_u_d[i]] = LOCAL_INDEX_GENERAL(_v[i], _u[i], nodes, symmetries);
    adjacency[_v[i]%based_nodes][_v_d[i]] = LOCAL_INDEX_GENERAL(_u[i], _v[i], nodes, symmetries);
  }
}

static bool check_isolated_vertices_general(const int u[2], const int v[2], const int nodes, const int symmetries,
					    const int degree, const int *num_degrees, int (*adjacency)[degree])
{
  int based_nodes = nodes/symmetries;
  if(check_isolated_vertex_general(u[0], based_nodes, degree, num_degrees, adjacency) ||
     check_isolated_vertex_general(v[0], based_nodes, degree, num_degrees, adjacency) ||
     check_isolated_vertex_general(u[1], based_nodes, degree, num_degrees, adjacency) ||
     check_isolated_vertex_general(v[1], based_nodes, degree, num_degrees, adjacency)){
    ODP_Restore_adjacency_general_s(nodes, degree, symmetries, adjacency);
    return false;
  }
  else{
    return true;
  }
}

void ODP_Restore_adjacency_general(const int nodes, const int degree, int (*adjacency)[degree])
{
  ODP_Restore_adjacency_general_s(nodes, degree, 1, adjacency);
}

static int get_degree_index_general(const int u, const int v, const int u_d, const int nodes,
				    const int symmetries, const int degree, const int *num_degrees,
				    const int (*adjacency)[degree])
{
  int based_nodes = nodes/symmetries;
  if(u == v){ // loop
    int d = (!num_degrees)? degree : num_degrees[v%based_nodes];
    for(int i=0;i<d;i++)
      if(GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, adjacency, v, i) == u && i != u_d)
	return i;
  }
  else if(symmetries%2 == 0 && abs(u-v) == nodes/2){
    return u_d;
  }
  else{
    int d = (!num_degrees)? degree : num_degrees[v%based_nodes];
    for(int i=0;i<d;i++)
      if(GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, adjacency, v, i) == u)
	return i;
  }

  ERROR("Something Wrong ! [id=4]\n");
  return -1; // dummy
}

static int get_degree_index_grid(const int u, const int v, const int u_d, const int width, const int height,
				 const int symmetries, const int degree, const int *num_degrees,
				 const int (*adjacency)[degree])
{
  if(u == v){ // loop
    int d = (!num_degrees)? degree : num_degrees[global2local_vertex_grid(v, width, height, symmetries)];
    for(int i=0;i<d;i++)
      if(ODP_GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, v, i) == u && i != u_d)
	return i;
  }
  else if(symmetries != 1 && u == ODP_ROTATE(v, width, height, symmetries, 180)){
    return u_d;
  }
  else{
    int d = (!num_degrees)? degree : num_degrees[global2local_vertex_grid(v, width, height, symmetries)];
    for(int i=0;i<d;i++){
      if(ODP_GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, v, i) == u)
	return i;
    }
  }
  
  ERROR("Something Wrong ! [id=5]\n");
  return -1; // dummy
}

static void backup_restore_adjacency(const int u[2], const int u_d[2], const int v[2], const int v_d[2])
{
  for(int i=0;i<2;i++){
    _u[i]   = u[i];
    _v[i]   = v[i];
    _u_d[i] = u_d[i];
    _v_d[i] = v_d[i];
  }
}

static bool mutate_adjacency_1opt_general_s(const int u, const int u_d, const int nodes, const int degree,
					    const int *num_degrees, const int symmetries,
					    int adjacency[nodes/symmetries][degree])
{
  int based_nodes = nodes/symmetries;
  int v = GLOBAL_ADJ_GENERAL(nodes, degree, symmetries, adjacency, u, u_d);
  if(IS_DIAMETER_GENERAL(u, v, nodes, symmetries)) return false;
  int v_d = get_degree_index_general(u, v, u_d, nodes, symmetries, degree, num_degrees, adjacency);
  int new_v;
  if((u-v)%based_nodes == 0){
    if(symmetries <= 4) return false;
    while(1){
      new_v = v + based_nodes * get_random(symmetries);
      if(u == new_v || v == new_v) continue;
      else if(NORM(u-v, nodes) == NORM(new_v-u, nodes)) continue;
      else if(nodes%2 == 1) break;
      else if(/* nodes%2 == 0 && */ (u-new_v)%(nodes/2) != 0) break;
    }
  }
  else{
    int rnd   = (symmetries%2 == 1)? get_random(symmetries-1) : get_random(symmetries);
    new_v = (rnd != symmetries-1)? v + based_nodes*(rnd+1) : u + based_nodes*(symmetries/2);
    //  int rnd   = get_random(symmetries-1);
    //  new_v = v + based_nodes*(rnd+1);
  }
  int new_u = v - (new_v - u);
  int tmp[2] = {LOCAL_INDEX_GENERAL(new_v, u, nodes, symmetries), LOCAL_INDEX_GENERAL(new_u, v, nodes, symmetries)};
   if(!check_multiple_edges_general_s(u%based_nodes, u_d, tmp[0], nodes, degree, num_degrees, symmetries, adjacency) ||
      !check_multiple_edges_general_s(v%based_nodes, v_d, tmp[1], nodes, degree, num_degrees, symmetries, adjacency))
     return false;
   
  adjacency[u%based_nodes][u_d] = tmp[0];
  adjacency[v%based_nodes][v_d] = tmp[1];
  return true;
}

static bool check_rotated_edges_overlap_general(const int u0, const int v0, const int u1, const int v1,
                                                const int nodes, const int symmetries)
{
  int based_nodes = nodes/symmetries;
  int diff0 = (u0 > v0)? v0 - u0 + nodes : v0 - u0;
  int diff1 = (u1 > v1)? v1 - u1 + nodes : v1 - u1;
  int diff2 = (u1 < v1)? u1 - v1 + nodes : u1 - v1;

  if(diff0 == diff1 && u0%based_nodes == u1%based_nodes)
    return true;
  else if(diff0 == diff2 && u0%based_nodes == v1%based_nodes)
    return true;
  else
    return false;
}

static bool mutate_adjacency_2opt_general_s(const int nodes, const int degree, const int *num_degrees,
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

  v_d[0] = get_degree_index_general(u[0], v[0], u_d[0], nodes, symmetries, degree, num_degrees, adjacency);
  v_d[1] = get_degree_index_general(u[1], v[1], u_d[1], nodes, symmetries, degree, num_degrees, adjacency);
  backup_restore_adjacency(u, u_d, v, v_d);
  
  if(IS_DIAMETER_GENERAL(u[0], v[0], nodes, symmetries) && IS_DIAMETER_GENERAL(u[1], v[1], nodes, symmetries)){
    if((u[0] - u[1])%based_nodes == 0)
      return false;

    int tmp[2];
    if(get_random(2)){
      tmp[0] = LOCAL_INDEX_GENERAL(u[1], u[0], nodes, symmetries);
      tmp[1] = LOCAL_INDEX_GENERAL(u[0], u[1], nodes, symmetries);
    }
    else{
      tmp[0] = LOCAL_INDEX_GENERAL(v[1], u[0], nodes, symmetries);
      tmp[1] = LOCAL_INDEX_GENERAL(v[0], u[1], nodes, symmetries);
    }
    if(!check_multiple_edges_general_s(u[0]%based_nodes, u_d[0], tmp[0], nodes, degree, num_degrees, symmetries, adjacency) || 
       !check_multiple_edges_general_s(u[1]%based_nodes, u_d[1], tmp[1], nodes, degree, num_degrees, symmetries, adjacency))
      return false;
    
    adjacency[u[0]%based_nodes][u_d[0]] = tmp[0];
    adjacency[u[1]%based_nodes][u_d[1]] = tmp[1];
    return true;
  }
  else if(IS_DIAMETER_GENERAL(u[0], v[0], nodes, symmetries) || IS_DIAMETER_GENERAL(u[1], v[1], nodes, symmetries)){
    if(IS_DIAMETER_GENERAL(u[1], v[1], nodes, symmetries)){
      SWAP(&u[0], &u[1]); SWAP(&u_d[0], &u_d[1]);
      SWAP(&v[0], &v[1]); SWAP(&v_d[0], &v_d[1]);
    }

    int opposite = nodes/2, rnd = get_random(4), tmp[4];
    if(rnd == 0){ // u[0]--v[1], u[1]--u[1]', v[0]--v[1]'
      tmp[0] = LOCAL_INDEX_GENERAL(v[1],          u[0], nodes, symmetries);
      tmp[1] = LOCAL_INDEX_GENERAL(u[1]+opposite, u[1], nodes, symmetries);
      tmp[2] = LOCAL_INDEX_GENERAL(v[1]+opposite, v[0], nodes, symmetries);
      tmp[3] = LOCAL_INDEX_GENERAL(u[0],          v[1], nodes, symmetries);
    }
    else if(rnd == 1){ // u[0]--v[1]', v[0]--v[1], u[1]--u[1]'
      tmp[0] = LOCAL_INDEX_GENERAL(v[1]+opposite, u[0], nodes, symmetries);
      tmp[1] = LOCAL_INDEX_GENERAL(u[1]+opposite, u[1], nodes, symmetries);
      tmp[2] = LOCAL_INDEX_GENERAL(v[1],          v[0], nodes, symmetries);
      tmp[3] = LOCAL_INDEX_GENERAL(v[0],          v[1], nodes, symmetries);
    }
    else if(rnd == 2){ // u[0]--u[1], v[0]--u[1]', v[1]--v[1]'
      tmp[0] = LOCAL_INDEX_GENERAL(u[1],          u[0], nodes, symmetries);
      tmp[1] = LOCAL_INDEX_GENERAL(u[0],          u[1], nodes, symmetries);
      tmp[2] = LOCAL_INDEX_GENERAL(u[1]+opposite, v[0], nodes, symmetries);
      tmp[3] = LOCAL_INDEX_GENERAL(v[1]+opposite, v[1], nodes, symmetries);
    }
    else if(rnd == 3){ // u[0]--u[1]', u[1]--v[0], v[1]--v[1]'
      tmp[0] = LOCAL_INDEX_GENERAL(u[1]+opposite, u[0], nodes, symmetries);
      tmp[1] = LOCAL_INDEX_GENERAL(v[0],          u[1], nodes, symmetries);
      tmp[2] = LOCAL_INDEX_GENERAL(u[1],          v[0], nodes, symmetries);
      tmp[3] = LOCAL_INDEX_GENERAL(v[1]+opposite, v[1], nodes, symmetries);
    }
    
    if(!check_multiple_edges_general_s(u[0]%based_nodes, u_d[0], tmp[0], nodes, degree, num_degrees, symmetries, adjacency) ||
       !check_multiple_edges_general_s(u[1]%based_nodes, u_d[1], tmp[1], nodes, degree, num_degrees, symmetries, adjacency) ||
       !check_multiple_edges_general_s(v[0]%based_nodes, v_d[0], tmp[2], nodes, degree, num_degrees, symmetries, adjacency) ||
       !check_multiple_edges_general_s(v[1]%based_nodes, v_d[1], tmp[3], nodes, degree, num_degrees, symmetries, adjacency))
      return false;
	  
    adjacency[u[0]%based_nodes][u_d[0]] = tmp[0];
    adjacency[u[1]%based_nodes][u_d[1]] = tmp[1];
    adjacency[v[0]%based_nodes][v_d[0]] = tmp[2];
    adjacency[v[1]%based_nodes][v_d[1]] = tmp[3];
    return true;
  }

  // Two selected edges are symmetrical
  if(check_rotated_edges_overlap_general(u[0], v[0], u[1], v[1], nodes, symmetries))
    return mutate_adjacency_1opt_general_s(u[0], u_d[0], nodes, degree, num_degrees, symmetries, adjacency);
  
  int tmp[4];
  if(get_random(2)){ // u[0]--v[1], v[0]--u[1]
    if(IS_DIAMETER_GENERAL(u[0], v[1], nodes, symmetries) || IS_DIAMETER_GENERAL(v[0], u[1], nodes, symmetries))
      return false;
    else if(check_rotated_edges_overlap_general(u[0], v[1], v[0], u[1], nodes, symmetries))
      return false;
    tmp[0] = LOCAL_INDEX_GENERAL(v[1], u[0], nodes, symmetries);
    tmp[1] = LOCAL_INDEX_GENERAL(v[0], u[1], nodes, symmetries);
    tmp[2] = LOCAL_INDEX_GENERAL(u[1], v[0], nodes, symmetries);
    tmp[3] = LOCAL_INDEX_GENERAL(u[0], v[1], nodes, symmetries);
  }
  else{ // u[0]--u[1], v[0]--v[1]
    if(IS_DIAMETER_GENERAL(u[0], u[1], nodes, symmetries) || IS_DIAMETER_GENERAL(v[0], v[1], nodes, symmetries))
      return false;
    else if(check_rotated_edges_overlap_general(u[0], u[1], v[0], v[1], nodes, symmetries))
      return false;
    tmp[0] = LOCAL_INDEX_GENERAL(u[1], u[0], nodes, symmetries);
    tmp[1] = LOCAL_INDEX_GENERAL(u[0], u[1], nodes, symmetries);
    tmp[2] = LOCAL_INDEX_GENERAL(v[1], v[0], nodes, symmetries);
    tmp[3] = LOCAL_INDEX_GENERAL(v[0], v[1], nodes, symmetries);
  }

  if(!check_multiple_edges_general_s(u[0]%based_nodes, u_d[0], tmp[0], nodes, degree, num_degrees, symmetries, adjacency) ||
     !check_multiple_edges_general_s(u[1]%based_nodes, u_d[1], tmp[1], nodes, degree, num_degrees, symmetries, adjacency) ||
     !check_multiple_edges_general_s(v[0]%based_nodes, v_d[0], tmp[2], nodes, degree, num_degrees, symmetries, adjacency) ||
     !check_multiple_edges_general_s(v[1]%based_nodes, v_d[1], tmp[3], nodes, degree, num_degrees, symmetries, adjacency))
    return false;
  
  adjacency[u[0]%based_nodes][u_d[0]] = tmp[0];
  adjacency[u[1]%based_nodes][u_d[1]] = tmp[1];
  adjacency[v[0]%based_nodes][v_d[0]] = tmp[2];
  adjacency[v[1]%based_nodes][v_d[1]] = tmp[3];
  return true;
}

void ODP_Mutate_adjacency_general_s(const int nodes, const int degree, const int *num_degrees,
				    const int symmetries, int adjacency[nodes/symmetries][degree])
{
  CHECK_SYMMETRIES(nodes, symmetries);
  
  while(1){
    if(mutate_adjacency_2opt_general_s(nodes, degree, num_degrees, symmetries, adjacency))
      break;
  }
}

static bool check_rotated_edges_overlap_grid(const int u0, const int v0, const int u1, const int v1,
					     const int width, const int height, const int symmetries)
{
  if(symmetries == 2){
    return (u0 == ODP_ROTATE(u1, width, height, symmetries, 180) && v0 == ODP_ROTATE(v1, width, height, symmetries, 180) ||
	    u0 == ODP_ROTATE(v1, width, height, symmetries, 180) && v0 == ODP_ROTATE(u1, width, height, symmetries, 180));
  }
  else if(symmetries == 4){
    return (u0 == ODP_ROTATE(u1, width, height, symmetries,  90) && v0 == ODP_ROTATE(v1, width, height, symmetries,  90) ||
	    u0 == ODP_ROTATE(u1, width, height, symmetries, 180) && v0 == ODP_ROTATE(v1, width, height, symmetries, 180) ||
	    u0 == ODP_ROTATE(u1, width, height, symmetries, 270) && v0 == ODP_ROTATE(v1, width, height, symmetries, 270) ||
	    u0 == ODP_ROTATE(v1, width, height, symmetries,  90) && v0 == ODP_ROTATE(u1, width, height, symmetries,  90) ||
	    u0 == ODP_ROTATE(v1, width, height, symmetries, 180) && v0 == ODP_ROTATE(u1, width, height, symmetries, 180) ||
	    u0 == ODP_ROTATE(v1, width, height, symmetries, 270) && v0 == ODP_ROTATE(u1, width, height, symmetries, 270));
  }
  return false; 
}

void ODP_Mutate_adjacency_general(const int nodes, const int degree, const int *num_degrees,
				  int adjacency[nodes][degree])
{
  ODP_Mutate_adjacency_general_s(nodes, degree, num_degrees, 1, adjacency);
}

static bool mutate_adjacency_1opt_grid_s(const int u, const int u_d, const int width, const int height, const int degree,
					 const int *num_degrees, const int length, const int symmetries, int (*adjacency)[degree])
{
  int v = ODP_GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, u, u_d);
  if(IS_DIAMETER_GRID(u, v, width, height, symmetries)) return false;
  else if(symmetries == 2 && u == ODP_ROTATE(v, width, height, symmetries, 180)) return false;
  else if(symmetries == 4 && (u == ODP_ROTATE(v, width, height, symmetries, 90) || u == ODP_ROTATE(v, width, height, symmetries, 270))) return false;
  int v_d = get_degree_index_grid(u, v, u_d, width, height, symmetries, degree, num_degrees, adjacency);

  int rnd = get_random(symmetries), tmp[2];
  if(rnd != symmetries-1){
    int new_v = ODP_ROTATE(v, width, height, symmetries, (rnd+1)*(360/symmetries));
    int new_u = ODP_ROTATE(u, width, height, symmetries, (symmetries-rnd-1)*(360/symmetries));
    if(!CHECK_LENGTH(u,new_v,height,length)) return false;
    if(!check_multiple_edges_grid_s(new_v, v_d, u, width, height, degree, num_degrees, symmetries, adjacency) ||
       !check_multiple_edges_grid_s(new_u, u_d, v, width, height, degree, num_degrees, symmetries, adjacency))
      return false;
    set_adjacency(new_v, v_d, u, width, height, degree, symmetries, adjacency);
    set_adjacency(new_u, u_d, v, width, height, degree, symmetries, adjacency);
  }
  else{
    int new_v = ODP_ROTATE(u, width, height, symmetries, 180);
    int new_u = ODP_ROTATE(v, width, height, symmetries, 180);
    if(!CHECK_LENGTH(u,new_v,height,length) || !CHECK_LENGTH(v,new_u,height,length)) return false;
    if(!check_multiple_edges_grid_s(new_v, u_d, u, width, height, degree, num_degrees, symmetries, adjacency) ||
       !check_multiple_edges_grid_s(new_u, v_d, v, width, height, degree, num_degrees, symmetries, adjacency))
      return false;
    
    set_adjacency(new_v, u_d, u, width, height, degree, symmetries, adjacency);
    set_adjacency(new_u, v_d, v, width, height, degree, symmetries, adjacency);
  }
  return true;
}

static bool mutate_adjacency_2opt_grid_s(const int width, const int height, const int degree, const int *num_degrees,
					 const int length, const int symmetries, int (*adjacency)[degree])
{
  int u[2], v[2], u_d[2], v_d[2], nodes = width*height, based_nodes = nodes/symmetries;

  while(1){
    u[0] = get_random(nodes);
    u[1] = get_random(nodes);
    if(u[0] == u[1]) continue;

    u_d[0] = (!num_degrees)? get_random(degree) : get_random(num_degrees[global2local_vertex_grid(u[0],width,height,symmetries)]);
    v[0] = ODP_GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, u[0], u_d[0]);
    if(v[0] == u[1]) continue;

    u_d[1] = (!num_degrees)? get_random(degree) : get_random(num_degrees[global2local_vertex_grid(u[1],width,height,symmetries)]);
    v[1] = ODP_GLOBAL_ADJ_GRID(width, height, degree, symmetries, adjacency, u[1], u_d[1]);
    if(v[1] == u[0] || v[0] == v[1]) continue;
    break;
  }

  v_d[0] = get_degree_index_grid(u[0], v[0], u_d[0], width, height, symmetries, degree, num_degrees, adjacency);
  v_d[1] = get_degree_index_grid(u[1], v[1], u_d[1], width, height, symmetries, degree, num_degrees, adjacency);
  backup_restore_adjacency(u, u_d, v, v_d);

  if(IS_DIAMETER_GRID(u[0], v[0], width, height, symmetries) && IS_DIAMETER_GRID(u[1], v[1], width, height, symmetries)){
    //    if(symmetries == 2 && u[0] == ODP_ROTATE(u[1], width, height, symmetries, 180)){
    //      return false; // This case is u[0] == v[1]
    //    }
    if(symmetries == 4 && (u[0] == ODP_ROTATE(u[1], width, height, symmetries, 90) || (u[0] == ODP_ROTATE(u[1], width, height, symmetries, 270))))
      return false;

    int tmp[2];
    if(get_random(2)){ // u[0]--u[1], v[0]--v[1]
      if(!CHECK_LENGTH(u[0], u[1], height, length)) return false;
      tmp[0] = u[1];
      tmp[1] = u[0];
    }
    else{ // u[0]--v[1], u[1]--v[0]
      if(!CHECK_LENGTH(u[0], v[1], height, length)) return false;
      tmp[0] = v[1];
      tmp[1] = v[0];
    }
    
    if(!check_multiple_edges_grid_s(u[0], u_d[0], tmp[0], width, height, degree, num_degrees, symmetries, adjacency) ||
       !check_multiple_edges_grid_s(u[1], u_d[1], tmp[1], width, height, degree, num_degrees, symmetries, adjacency))
      return false;
    
    set_adjacency(u[0], u_d[0], tmp[0], width, height, degree, symmetries, adjacency);
    set_adjacency(u[1], u_d[1], tmp[1], width, height, degree, symmetries, adjacency);
    return true;
  }
  else if(IS_DIAMETER_GRID(u[0], v[0], width, height, symmetries) || IS_DIAMETER_GRID(u[1], v[1], width, height, symmetries)){
    if(IS_DIAMETER_GRID(u[1], v[1], width, height, symmetries)){
      SWAP(&u[0], &u[1]); SWAP(&u_d[0], &u_d[1]);
      SWAP(&v[0], &v[1]); SWAP(&v_d[0], &v_d[1]);
    }

    int rnd = get_random(4), tmp[4];
    int u1_opposite = ODP_ROTATE(u[1], width, height, symmetries, 180);
    int v1_opposite = ODP_ROTATE(v[1], width, height, symmetries, 180);
    if(rnd == 0){ // u[0]--v[1], u[1]--u[1]', v[0]--v[1]'
      if(!CHECK_LENGTH(u[0], v[1], height, length) || !CHECK_LENGTH(u[1], u1_opposite, height, length))
	return false;
      tmp[0] = v[1];
      tmp[1] = u1_opposite;
      tmp[2] = v1_opposite;
      tmp[3] = u[0];
    }
    else if(rnd == 1){ // u[0]--v[1]', v[0]--v[1], u[1]--u[1]'
      if(!CHECK_LENGTH(u[0], v1_opposite, height, length) || !CHECK_LENGTH(u[1], u1_opposite, height, length))
        return false;
      tmp[0] = v1_opposite;
      tmp[1] = u1_opposite;
      tmp[2] = v[1];
      tmp[3] = v[0];
    }
    else if(rnd == 2){ // u[0]--u[1], v[0]--u[1]', v[1]--v[1]'
      if(!CHECK_LENGTH(u[0], u[1], height, length) || !CHECK_LENGTH(v[1], v1_opposite, height, length))
        return false;
      tmp[0] = u[1];
      tmp[1] = u[0];
      tmp[2] = u1_opposite;
      tmp[3] = v1_opposite;
    }
    else if(rnd == 3){ // u[0]--u[1]', u[1]--v[0], v[1]--v[1]'
      if(!CHECK_LENGTH(u[0], u1_opposite, height, length) || !CHECK_LENGTH(v[1], v1_opposite, height, length))
        return false;
      tmp[0] = u1_opposite;
      tmp[1] = v[0];
      tmp[2] = u[1];
      tmp[3] = v1_opposite;
    }

     if(!check_multiple_edges_grid_s(u[0], u_d[0], tmp[0], width, height, degree, num_degrees, symmetries, adjacency) ||
	!check_multiple_edges_grid_s(u[1], u_d[1], tmp[1], width, height, degree, num_degrees, symmetries, adjacency) ||
	!check_multiple_edges_grid_s(v[0], v_d[0], tmp[2], width, height, degree, num_degrees, symmetries, adjacency) ||
	!check_multiple_edges_grid_s(v[1], v_d[1], tmp[3], width, height, degree, num_degrees, symmetries, adjacency))
       return false;
     
    set_adjacency(u[0], u_d[0], tmp[0], width, height, degree, symmetries, adjacency);
    set_adjacency(u[1], u_d[1], tmp[1], width, height, degree, symmetries, adjacency);
    set_adjacency(v[0], v_d[0], tmp[2], width, height, degree, symmetries, adjacency);
    set_adjacency(v[1], v_d[1], tmp[3], width, height, degree, symmetries, adjacency);
    return true;
  }

   // Two selected edges are symmetrical
  if(check_rotated_edges_overlap_grid(u[0], v[0], u[1], v[1], width, height, symmetries))
    return mutate_adjacency_1opt_grid_s(u[0], u_d[0], width, height, degree, num_degrees, length, symmetries, adjacency);

  int tmp[4];
  if(get_random(2)){ // u[0]--v[1], v[0]--u[1]
    if(!CHECK_LENGTH(u[0], v[1], height, length) || !CHECK_LENGTH(v[0], u[1], height, length))
      return false;
    else if(IS_DIAMETER_GRID(u[0], v[1], width, height, symmetries) || IS_DIAMETER_GRID(v[0], u[1], width, height, symmetries))
      return false;
    else if(check_rotated_edges_overlap_grid(u[0], v[1], v[0], u[1], width, height, symmetries))
      return false;
    tmp[0] = v[1];
    tmp[1] = v[0];
    tmp[2] = u[1];
    tmp[3] = u[0];
  }
  else{ // u[0]--u[1], v[0]--v[1]
    if(!CHECK_LENGTH(u[0], u[1], height, length) || !CHECK_LENGTH(v[0], v[1], height, length))
      return false;
    else if(IS_DIAMETER_GRID(u[0], u[1], width, height, symmetries) || IS_DIAMETER_GRID(v[0], v[1], width, height, symmetries))
      return false;
    else if(check_rotated_edges_overlap_grid(u[0], u[1], v[0], v[1], width, height, symmetries))
      return false;
    tmp[0] = u[1];
    tmp[1] = u[0];
    tmp[2] = v[1];
    tmp[3] = v[0];
  }

  if(!check_multiple_edges_grid_s(u[0], u_d[0], tmp[0], width, height, degree, num_degrees, symmetries, adjacency) ||
     !check_multiple_edges_grid_s(u[1], u_d[1], tmp[1], width, height, degree, num_degrees, symmetries, adjacency) ||
     !check_multiple_edges_grid_s(v[0], v_d[0], tmp[2], width, height, degree, num_degrees, symmetries, adjacency) ||
     !check_multiple_edges_grid_s(v[1], v_d[1], tmp[3], width, height, degree, num_degrees, symmetries, adjacency))
    return false;
    
  set_adjacency(u[0], u_d[0], tmp[0], width, height, degree, symmetries, adjacency);
  set_adjacency(u[1], u_d[1], tmp[1], width, height, degree, symmetries, adjacency);
  set_adjacency(v[0], v_d[0], tmp[2], width, height, degree, symmetries, adjacency);
  set_adjacency(v[1], v_d[1], tmp[3], width, height, degree, symmetries, adjacency);
  return true;
}

void ODP_Mutate_adjacency_grid_s(const int width, const int height, const int degree,
				 const int *num_degrees, const int length, const int symmetries,
				 int (*adjacency)[degree])
{
  CHECK_PARAMETERS(width*height, degree);
  CHECK_SYMMETRIES_GRID(symmetries);  
  CHECK_SYMMETRIES_WH(symmetries, width, height);
  
  while(1){
    if(mutate_adjacency_2opt_grid_s(width, height, degree, num_degrees, length, symmetries, adjacency))
      break;
  }
}

void ODP_Conv_adjacency2edge_general_s(const int nodes, const int degree, const int *num_degrees,
                                       const int *adjacency, const int symmetries, int (*edge)[2])
{
  CHECK_SYMMETRIES(nodes, symmetries);

  int (*tmp)[degree] = malloc(sizeof(int) * nodes * degree);
  int based_nodes = nodes/symmetries;
  for(int i=0;i<symmetries;i++){
    for(int j=0;j<based_nodes;j++){
      int d = (!num_degrees)? degree : num_degrees[j];
      for(int k=0;k<d;k++){
	int v = *(adjacency + j * degree + k) + i * based_nodes;
        tmp[i*based_nodes+j][k] = NORM(v, nodes);
      }
    }
  }

  ODP_Conv_adjacency2edge(nodes, degree, num_degrees, (int *)tmp, edge);
  free(tmp);
}

void ODP_Conv_edge2adjacency_general_s(const int nodes, const int lines, const int degree, const int edge[lines][2],
                                       const int symmetries, int (*adjacency)[degree])
{
  CHECK_SYMMETRIES(nodes, symmetries);

  int based_nodes = nodes/symmetries;
  int num_degrees[based_nodes];
  for(int i=0;i<based_nodes;i++)
    num_degrees[i] = 0;

  for(int i=0;i<nodes/symmetries;i++)
    for(int j=0;j<degree;j++)
      adjacency[i][j] = NOT_DEFINED;

  for(int i=0;i<lines;i++){
    int n1 = edge[i][0];
    int n2 = edge[i][1];
    if(n1 < based_nodes)
      adjacency[n1][num_degrees[n1]++] = n2;
    if(n2 < based_nodes)
      adjacency[n2][num_degrees[n2]++] = n1;
  }
}

void ODP_Srand(const unsigned int seed)
{
  srand(seed);
}

void ODP_Generate_random_general(const int nodes, const int degree, int (*edge)[2])
{
  CHECK_PARAMETERS(nodes, degree);

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
  ODP_Conv_edge2adjacency(nodes, lines, degree, edge, adjacency);

  // Give randomness
  for(int i=0;i<lines*GEN_GRAPH_ITERS;i++)
    ODP_Mutate_adjacency_general(nodes, degree, NULL, (int (*)[degree])adjacency);

  // Repeat until there are no unreachable vertices
  while(simple_bfs(nodes, degree, adjacency))
    ODP_Mutate_adjacency_general(nodes, degree, NULL, (int (*)[degree])adjacency);

  ODP_Conv_adjacency2edge(nodes, degree, NULL, adjacency, edge);
  free(adjacency);
}

void ODP_Generate_random_general_s(const int nodes, const int degree, const int symmetries, int (*edge)[2])
{
  CHECK_SYMMETRIES(nodes, symmetries);
  CHECK_PARAMETERS(nodes, degree);

  if(symmetries == 1){
    ODP_Generate_random_general(nodes, degree, edge);
    return;
  }
  
  int based_nodes = nodes/symmetries;
  int lines       = (nodes*degree)/2;
  int based_lines = lines/symmetries;

  if(based_nodes%2==0 || degree%2 == 0){
    ODP_Generate_random_general(based_nodes, degree, edge);
    for(int i=1;i<symmetries;i++){
      for(int j=0;j<based_lines;j++){
	for(int k=0;k<2;k++){
	  int v = edge[j][k] + based_nodes * i;
	  edge[i*based_lines+j][k] = NORM(v, nodes);
	}
      }
    }
  }
  else{
    ODP_Generate_random_general(based_nodes, degree-1, edge);
    int based_lines_shrink = (based_nodes*(degree-1))/2;
    for(int i=1;i<symmetries;i++){
      for(int j=0;j<based_lines_shrink;j++){
	for(int k=0;k<2;k++){
	  int v = edge[j][k] + based_nodes * i;
	  edge[i*based_lines_shrink+j][k] = NORM(v, nodes);
	}
      }
    }
    int offset = lines - nodes/2;
    for(int i=0;i<nodes/2;i++){
      edge[offset+i][0] = i;
      edge[offset+i][1] = i + nodes/2;
    }
  }

  int *adjacency = malloc(sizeof(int) * based_nodes * degree);
  ODP_Conv_edge2adjacency_general_s(nodes, lines, degree, edge, symmetries, (int (*)[degree])adjacency);

  // Give randomness
  for(int i=0;i<based_lines*GEN_GRAPH_ITERS;i++)
    ODP_Mutate_adjacency_general_s(nodes, degree, NULL, symmetries, (int (*)[degree])adjacency);

  // Repeat until there are no unreachable vertices
  // (omit)

  ODP_Conv_adjacency2edge_general_s(nodes, degree, NULL, adjacency, symmetries, edge);
  free(adjacency);
}

void ODP_Conv_edge2adjacency_grid_s(const int width, const int height, const int lines, const int degree,
				    const int edge[lines][2], const int symmetries, int (*adjacency)[degree])
{
  CHECK_SYMMETRIES_GRID(symmetries);
  CHECK_SYMMETRIES_WH(symmetries, width, height);
  
  int nodes = width * height;
  if(symmetries == 1){
    ODP_Conv_edge2adjacency(nodes, lines, degree, edge, (int *)adjacency);
    return;
  }

  int based_nodes = nodes/symmetries;
  int num_degrees[based_nodes];
  for(int i=0;i<based_nodes;i++)
    num_degrees[i] = 0;
  
  for(int i=0;i<nodes/symmetries;i++)
    for(int j=0;j<degree;j++)
      adjacency[i][j] = NOT_DEFINED;

  int based_width = width/2;
  int based_height = (symmetries == 2)? height : height/2;
  for(int i=0;i<lines;i++){
    int n1 = edge[i][0];
    int n2 = edge[i][1];    
    if(WIDTH(n1,height)<based_width && HEIGHT(n1,height)<based_height){
      int local_n1 = ODP_LOCAL_INDEX_GRID(n1, width, height, symmetries);
      adjacency[local_n1][num_degrees[local_n1]++] = n2;
    }
    if(WIDTH(n2,height)<based_width && HEIGHT(n2,height)<based_height){
      int local_n2 = ODP_LOCAL_INDEX_GRID(n2, width, height, symmetries);
      adjacency[local_n2][num_degrees[local_n2]++] = n1;
    }
  }
}

void ODP_Write_edge_general(const int lines, const int edge[lines][2], char *fname)
{
  FILE *fp = NULL;

  if((fp = fopen(fname, "w")) == NULL)
    ERROR("Cannot open %s\n", fname);

  for(int i=0;i<lines;i++)
    fprintf(fp, "%d %d\n", edge[i][0], edge[i][1]);

  fclose(fp);
}

void ODP_Write_edge_grid(const int lines, const int height, const int edge[lines][2], char *fname)
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

void ODP_Malloc(uint64_t **a, const size_t s, const bool enable_avx2)
{
  if(enable_avx2)
    *a = _mm_malloc(s, ALIGN_VALUE);
  else
    posix_memalign((void **)a, ALIGN_VALUE, s);
}

void ODP_Free(uint64_t *a, const bool enable_avx2)
{
  if(enable_avx2)
    _mm_free(a);
  else
    free(a);
}

double ODP_Get_mem_usage(const int kind, const int nodes, const int degree, const int symmetries,
                         const int *num_degrees, const int procs, const bool is_cpu)
{
  int Mbyte = 1024*1024;
  int chunk = (is_cpu)? CPU_CHUNK : GPU_CHUNK;
  double AB_mem = (kind == ASPL_NORMAL)? (nodes*((double)nodes/(4*symmetries*procs))) : (double)16*nodes*chunk;

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

int ODP_Get_kind(const int nodes, const int degree, const int* num_degrees, const int symmetries,
                 const int procs, const int is_cpu)
{
  int kind;
  char *val = getenv("ODP_ASPL");
  if(!val){
    double normal_mem_usage = ODP_Get_mem_usage(ASPL_NORMAL, nodes, degree, symmetries, num_degrees, procs, is_cpu);
    if(normal_mem_usage <= MEM_THRESHOLD)
      kind = ASPL_NORMAL;
    else
      kind = ASPL_SAVING;
  }
  else if(strcmp(val, "NORMAL") == 0){
    kind = ASPL_NORMAL;
  }
  else if(strcmp(val, "SAVING") == 0){
    kind = ASPL_SAVING;
  }
  else{
    ERROR("Unknown ASPL value (%s)\n", val);
  }

  return kind;
}

void ODP_Profile(const char* name, const int kind, const int symmetries, const double mem_usage,
                 const double elapsed_time, const unsigned int times, const int procs)
{
  char kind_name[7], hostname[MAX_HOSTNAME_LENGTH];
  if(kind == ASPL_NORMAL) strcpy(kind_name, "NORMAL");
  else                    strcpy(kind_name, "SAVING");
  gethostname(hostname, sizeof(hostname));
  time_t t = time(NULL);

  printf("------ Profile for SET_ASPL ------\n");
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

bool ODP_Check_loop(const int lines, int edge[lines][2])
{
  for(int i=0;i<lines;i++)
    if(edge[i][0] == edge[i][1]){
      printf("Loop is found in line %d\n", i+1);
      return true;
    }

  return false;
}

bool ODP_Check_multiple_edges(const int lines, int edge[lines][2])
{
  for(int i=0;i<lines;i++)
    for(int j=i+1;j<lines;j++)
      if(has_multiple_edges(edge[i][0], edge[i][1], edge[j][0], edge[j][1])){
        printf("Multiple edeges are found in lines %d %d\n", i+1, j+1);
        return true;
      }

  return false;
}

int ODP_Get_length(const int lines, const int edge[lines][2], const int height)
{
  int length = 0;
  for(int i=0;i<lines;i++)
    length = MAX(length, DISTANCE(edge[i][0], edge[i][1], height));

  return length;
}

bool ODP_Check_general(char *fname)
{
  FILE *fp;
  if((fp = fopen(fname, "r")) == NULL)
    ERROR("File not found\n");

  int n1=1, n2=-1;
  fscanf(fp, "%d %d", &n1, &n2);
  fclose(fp);

  return (n2 != -1)? true : false;
}

int ODP_Get_degree(const int nodes, const int lines, const int edge[lines][2])
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

int ODP_Get_nodes(const int lines, const int (*edge)[2])
{
  int max = 0;
  for(int i=0;i<lines;i++){
    max = MAX(max, edge[i][0]);
    max = MAX(max, edge[i][1]);
  }

  return max + 1;
}

int ODP_Get_lines(const char* fname)
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

void ODP_Read_edge_general(const char* fname, int (*edge)[2])
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

void ODP_Read_edge_grid(const char *fname, int *w, int *h, int (*edge)[2])
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
void ODP_Set_lbounds_general(const int nodes, const int degree, int *low_diameter, double *low_ASPL)
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

// This function is inherited from "http://research.nii.ac.jp/graphgolf/pl/lower-lattice.pl".
void ODP_Set_lbounds_grid(const int m, const int n, const int degree, const int length, int *low_diameter, double *low_ASPL)
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

void ODP_Conv_adjacency2edge(const int nodes, const int degree, const int *num_degrees,
                             const int *adjacency, int (*edge)[2])
{
  char (*tmp)[degree] = malloc(sizeof(char) * nodes * degree);
  for(int i=0;i<nodes;i++)
    for(int j=0;j<degree;j++)
      tmp[i][j] = NOT_VISITED;

  int lines = get_lines(nodes, degree, num_degrees);
  int j = 0;
  for(int u=0;u<nodes;u++){
    int d = (!num_degrees)? degree : num_degrees[u];
    for(int i=0;i<d;i++){
      int v = *(adjacency + u * degree + i);
      if(tmp[u][i] == NOT_VISITED){
        if(j >= lines) ERROR("Something Wrong ! [id=0]\n");
        edge[j][0] = u;
        edge[j][1] = v;
        tmp[u][i]  = VISITED;
        j++;
        int k=0;
	int d2 = (!num_degrees)? degree : num_degrees[v];
        for(k=0;k<d2;k++){
          if(*(adjacency + v * degree + k) == u && tmp[v][k] == NOT_VISITED){
            tmp[v][k] = VISITED;
            break;
          }
        }
        if(k==d2)
          ERROR("Something Wrong ! [id=1]\n");
      }
    }
  }

  free(tmp);
}

void ODP_Conv_edge2adjacency(const int nodes, const int lines, const int degree, const int edge[lines][2],
                             int *adjacency) // int adjacency[nodes][degree]
{
  ODP_Conv_edge2adjacency_general_s(nodes, lines, degree, edge, 1, (int (*)[degree])adjacency);
}

void ODP_Set_degrees(const int nodes, const int lines, int edge[lines][2],
                     int* num_degrees)
{
  for(int i=0;i<nodes;i++)
    num_degrees[i] = 0;

  for(int i=0;i<lines;i++){
    num_degrees[edge[i][0]]++;
    num_degrees[edge[i][1]]++;
  }
}

void ODP_Mutate_adjacency_grid(const int width, const int height, const int degree,
			       const int *num_degrees, const int length,
			       int (*adjacency)[degree])
{
  CHECK_PARAMETERS(width*height, degree);
  while(1){
    if(mutate_adjacency_2opt_grid_s(width, height, degree, num_degrees, length, 1, adjacency))
      break;
  }
}

// Inherited from http://research.nii.ac.jp/graphgolf/c/create-lattice.c
void ODP_Generate_random_grid(const int width, const int height, const int degree, const int length,
                              int (*edge)[2])
{
  int nodes = width * height;
  CHECK_PARAMETERS(nodes, degree);

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
  int *adjacency = malloc(sizeof(int)*nodes*degree);
  ODP_Conv_edge2adjacency(nodes, lines, degree, edge, adjacency);

  // Give randomness
  for(int i=0;i<lines*GEN_GRAPH_ITERS;i++)
    ODP_Mutate_adjacency_grid(width, height, degree, NULL, length, (int (*)[degree])adjacency);

  // Repeat until there are no unreachable vertices
  while(simple_bfs(nodes, degree, adjacency))
    ODP_Mutate_adjacency_grid(width, height, degree, NULL, length, (int (*)[degree])adjacency);

  ODP_Conv_adjacency2edge(nodes, degree, NULL, adjacency, edge);
  free(adjacency);
}

void ODP_Generate_random_grid_s(const int width, const int height, const int degree, const int length,
                                const int symmetries, int (*edge)[2])
{
  CHECK_PARAMETERS(width*height, degree);
  CHECK_SYMMETRIES_GRID(symmetries);
  CHECK_SYMMETRIES_WH(symmetries, width, height);

  if(symmetries == 1){
    ODP_Generate_random_grid(width, height, degree, length, edge);
    return;
  }

  int nodes        = width * height;
  int lines        = (nodes * degree) / 2;
  int based_lines  = lines / symmetries;
  int based_nodes  = nodes / symmetries;
  int based_width  = width / 2;
  int based_height = (symmetries == 2)? height : height/2;

  if(based_nodes%2==0 || degree%2 == 0){
    ODP_Generate_random_grid(based_width, based_height, degree, length, edge);

    for(int i=0;i<based_lines;i++)
      for(int j=0;j<2;j++)
        edge[i][j] = WIDTH(edge[i][j], based_height) * height + HEIGHT(edge[i][j], based_height);

    if(symmetries == 2){
      for(int i=0;i<based_lines;i++)
        for(int j=0;j<2;j++)
          edge[based_lines+i][j] = ODP_ROTATE(edge[i][j], width, height, symmetries, 180);
    }
    else{ // symmetries == 4
      for(int i=0;i<based_lines;i++){
        for(int j=0;j<2;j++){
          edge[based_lines  +i][j] = ODP_ROTATE(edge[i][j], width, height, symmetries,  90);
          edge[based_lines*2+i][j] = ODP_ROTATE(edge[i][j], width, height, symmetries, 180);
          edge[based_lines*3+i][j] = ODP_ROTATE(edge[i][j], width, height, symmetries, 270);
        }
      }
    }
  }
  else{
    ODP_Generate_random_grid(based_width, based_height, degree-1, length, edge);
    int based_lines_shrink = (based_nodes*(degree-1))/2;

    for(int i=0;i<based_lines_shrink;i++)
      for(int j=0;j<2;j++)
        edge[i][j] = WIDTH(edge[i][j], based_height) * height + HEIGHT(edge[i][j], based_height);

    if(symmetries == 2){
      for(int i=0;i<based_lines_shrink;i++)
        for(int j=0;j<2;j++)
          edge[based_lines_shrink+i][j] = ODP_ROTATE(edge[i][j], width, height, symmetries, 180);

      int k = based_lines_shrink*symmetries;
      for(int w=0;w<width;w+=2){
        for(int h=0;h<height;h++){
          edge[k][0] = h + w * height;
          edge[k][1] = edge[k][0] + height;
          k++;
        }
      }
    }
    else{ // symmetries == 4
      for(int i=0;i<based_lines_shrink;i++){
        for(int j=0;j<2;j++){
          edge[based_lines_shrink  +i][j] = ODP_ROTATE(edge[i][j], width, height, symmetries,  90);
          edge[based_lines_shrink*2+i][j] = ODP_ROTATE(edge[i][j], width, height, symmetries, 180);
          edge[based_lines_shrink*3+i][j] = ODP_ROTATE(edge[i][j], width, height, symmetries, 270);
        }
      }
      int k = based_lines_shrink*symmetries;
      for(int w=0;w<based_width/2;w++){
        for(int h=0;h<based_height;h++){
          edge[k  ][0] = w * 2 * height + h;
          edge[k  ][1] = edge[k][0] + height;
          edge[k+1][0] = ODP_ROTATE(edge[k][0], width, height, symmetries,  90);
          edge[k+1][1] = ODP_ROTATE(edge[k][1], width, height, symmetries,  90);
          edge[k+2][0] = ODP_ROTATE(edge[k][0], width, height, symmetries, 180);
          edge[k+2][1] = ODP_ROTATE(edge[k][1], width, height, symmetries, 180);
          edge[k+3][0] = ODP_ROTATE(edge[k][0], width, height, symmetries, 270);
          edge[k+3][1] = ODP_ROTATE(edge[k][1], width, height, symmetries, 270);
          k+=4;
        }
      }
      if(based_width%2==1){
        for(int h=0;h<based_height/2;h++){
          edge[k  ][0] = (based_width - 1) * height + h * 2;
          edge[k  ][1] = edge[k][0] + 1;
          edge[k+1][0] = ODP_ROTATE(edge[k][0], width, height, symmetries,  90);
          edge[k+1][1] = ODP_ROTATE(edge[k][1], width, height, symmetries,  90);
          edge[k+2][0] = ODP_ROTATE(edge[k][0], width, height, symmetries, 180);
          edge[k+2][1] = ODP_ROTATE(edge[k][1], width, height, symmetries, 180);
          edge[k+3][0] = ODP_ROTATE(edge[k][0], width, height, symmetries, 270);
          edge[k+3][1] = ODP_ROTATE(edge[k][1], width, height, symmetries, 270);
          k+=4;
        }
        if(based_height%2==1){
          if(length <= 1)
            ERROR("Graphs with odd width/2 and length/2 and length 1 cannot be created.");
          edge[k  ][0] = (based_width - 1) * height + (based_height - 1);
          edge[k  ][1] = based_width * height + (based_height);
          edge[k+1][0] = (based_width - 1) * height + based_height;
          edge[k+1][1] = based_width * height + (based_height - 1);
        }
      }
    }
  }
  int (*adjacency)[degree] = malloc(sizeof(int) * based_nodes * degree);
  ODP_Conv_edge2adjacency_grid_s(width, height, lines, degree, edge, symmetries, adjacency);

  // Give randomness
  for(int i=0;i<based_lines*GEN_GRAPH_ITERS;i++)
    ODP_Mutate_adjacency_grid_s(width, height, degree, NULL, length, symmetries, adjacency);

  // Repeat until there are no unreachable vertices
  // (omit)

  ODP_Conv_adjacency2edge_grid_s(width, height, degree, NULL, adjacency, symmetries, edge);

  free(adjacency);
}

void ODP_Create_itable(const int width, const int height, const int symmetries, int *itable)
{
  int based_nodes = (width*height)/symmetries;
  if(symmetries == 2){
    for(int i=0;i<based_nodes;i++){
      itable[i] = i;
      itable[ODP_ROTATE(i, width, height, symmetries, 180)] = i + based_nodes;
    }
  }
  else if(symmetries == 4){
    int based_height = height/2;
    for(int i=0;i<based_nodes;i++){
      int v = (i/based_height) * height + (i%based_height);
      itable[v] = i;
      itable[ODP_ROTATE(v, width, height, symmetries,  90)] = i + based_nodes;
      itable[ODP_ROTATE(v, width, height, symmetries, 180)] = i + based_nodes * 2;
      itable[ODP_ROTATE(v, width, height, symmetries, 270)] = i + based_nodes * 3;
    }
  }
  else
    ERROR("Something Wrong !");
}

