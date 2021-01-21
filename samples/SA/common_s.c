#include "common_s.h"
static int _u[2], _v[2], _u_d[2], _v_d[2], _nodes, _symmetries, _kind;

double get_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1.0e-6 * t.tv_usec;
}

static int NORM(int x, const int nodes)
{
  while(x < 0 || x >= nodes)
    x = (x < 0)? x + nodes : x - nodes;

  return x;
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

// return adjacency[global_vertex][d];
static int GLOBAL_ADJ(const int nodes, const int degree, const int symmetries,
                      const int (*adjacency)[degree], const int global_vertex, const int d)
{
  int based_nodes = nodes/symmetries;
  int n = adjacency[global_vertex%based_nodes][d] + (global_vertex/based_nodes)*based_nodes;
  return NORM(n, nodes);
}

static int LOCAL_VERTEX(const int global_vertex, const int position, const int nodes, const int symmetries)
{
  int based_nodes = nodes/symmetries;
  return NORM(global_vertex - (position/based_nodes)*based_nodes, nodes);
}

static void backup_restore_adjacency(const int u[2], const int u_d[2], const int v[2], const int v_d[2],
				     const int nodes, const int symmetries, const int kind)
{
  int based_nodes = nodes/symmetries;
  memcpy(  _u,   u, sizeof(int)*2);
  memcpy(  _v,   v, sizeof(int)*2);
  memcpy(_u_d, u_d, sizeof(int)*2);
  memcpy(_v_d, v_d, sizeof(int)*2);
  _nodes      = nodes;
  _symmetries = symmetries;
  _kind       = kind;
}

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

bool accept_s(const int nodes, const int current_diameter, const int diameter,
              const double current_ASPL, const double ASPL, const double temp, const int symmetries)
{
  if(diameter < current_diameter){
    return true;
  }
  else if(diameter > current_diameter){
    return false;
  }
  else{ //  diameter == current_diameter
    if(ASPL <= current_ASPL){
      return true;
    }
    else{
      double diff = ((current_ASPL-ASPL)*nodes*(nodes-1))/symmetries;
      if(exp(diff/temp) > uniform_rand()){
        return true;
      }
      else{
        return false;
      }
    }
  }
}

static int get_random(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
}

void restore_adjacency(const int degree, int *adjacency)
{
  int based_nodes = _nodes/_symmetries;
  adjacency[_u[0]%based_nodes * degree + _u_d[0]] = LOCAL_VERTEX(_v[0], _u[0], _nodes, _symmetries);
  adjacency[_v[0]%based_nodes * degree + _v_d[0]] = LOCAL_VERTEX(_u[0], _v[0], _nodes, _symmetries);
  if(_kind == MUTATE_1OPT) return;
  adjacency[_u[1]%based_nodes * degree + _u_d[1]] = LOCAL_VERTEX(_v[1], _u[1], _nodes, _symmetries);
  adjacency[_v[1]%based_nodes * degree + _v_d[1]] = LOCAL_VERTEX(_u[1], _v[1], _nodes, _symmetries);
}

static int get_degree_index(const int u, const int v, const int u_d, const int nodes,
			    const int degree,  const int symmetries, const int (*adjacency)[degree])
{
  int based_nodes = nodes/symmetries;
  if(u == v){ // loop
    for(int i=0;i<degree;i++)
      if(GLOBAL_ADJ(nodes, degree, symmetries, adjacency, v, i) == u && i != u_d)
	return i;
  }
  else if(symmetries%2 == 0 && abs(u-v) == nodes/2){
    return u_d;
  }
  else{
    for(int i=0;i<degree;i++)
      if(GLOBAL_ADJ(nodes, degree, symmetries, adjacency, v, i) == u)
        return i;
  }

  ERROR("Something wrong ! [id=222]\n");
  return -1; // dummy
}

static bool check_isolated_vertex(const int x, const int degree, const int (*adjacency)[degree])
{
  for(int i=1;i<degree;i++)
    if(adjacency[x][0] != adjacency[x][i])
      return false;

  return true;
}

static bool mutate_adjacency_2opt_s(const int nodes, const int degree,
                                    const int symmetries, int adjacency[nodes/symmetries][degree])
{
  int u[2], v[2], u_d[2], v_d[2], based_nodes = nodes/symmetries;

  while(1){
    u[0] = get_random(nodes);
    u[1] = get_random(nodes);
    if(u[0] == u[1]) continue;

    u_d[0] = get_random(degree);
    v[0] = GLOBAL_ADJ(nodes, degree, symmetries, adjacency, u[0], u_d[0]);
    if(v[0] == u[1]) continue;

    u_d[1] = get_random(degree);
    v[1] = GLOBAL_ADJ(nodes, degree, symmetries, adjacency, u[1], u_d[1]);
    if(v[1] == u[0] || v[0] == v[1]) continue;
    break;
  }

  v_d[0] = get_degree_index(u[0], v[0], u_d[0], nodes, degree, symmetries, adjacency);
  v_d[1] = get_degree_index(u[1], v[1], u_d[1], nodes, degree, symmetries, adjacency);
  backup_restore_adjacency(u, u_d, v, v_d, nodes, symmetries, MUTATE_2OPT);

  if(IS_DIAMETER(u[0], v[0], nodes, symmetries) && IS_DIAMETER(u[1], v[1], nodes, symmetries)){
    if((u[0] - u[1])%based_nodes == 0){
      return false;
    }
    else{
      adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_VERTEX(u[1], u[0], nodes, symmetries);
      adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_VERTEX(u[0], u[1], nodes, symmetries);
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
      adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_VERTEX(v[1],          u[0], nodes, symmetries);
      adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_VERTEX(u[1]+opposite, u[1], nodes, symmetries);
      adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_VERTEX(v[1]+opposite, v[0], nodes, symmetries);
      adjacency[v[1]%based_nodes][v_d[1]] = LOCAL_VERTEX(u[0],          v[1], nodes, symmetries);
    }
    else if(rnd == 1){
      adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_VERTEX(v[1]+opposite, u[0], nodes, symmetries);
      adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_VERTEX(u[1]+opposite, u[1], nodes, symmetries);
      adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_VERTEX(v[1],          v[0], nodes, symmetries);
      adjacency[v[1]%based_nodes][v_d[1]] = LOCAL_VERTEX(v[0],          v[1], nodes, symmetries);
    }
    else if(rnd == 2){
      adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_VERTEX(u[1],          u[0], nodes, symmetries);
      adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_VERTEX(u[0],          u[1], nodes, symmetries);
      adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_VERTEX(u[1]+opposite, v[0], nodes, symmetries);
      adjacency[v[1]%based_nodes][v_d[1]] = LOCAL_VERTEX(v[1]+opposite, v[1], nodes, symmetries);
    }
    else if(rnd == 3){
      adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_VERTEX(u[1]+opposite, u[0], nodes, symmetries);
      adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_VERTEX(v[0],          u[1], nodes, symmetries);
      adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_VERTEX(u[1],          v[0], nodes, symmetries);
      adjacency[v[1]%based_nodes][v_d[1]] = LOCAL_VERTEX(v[1]+opposite, v[1], nodes, symmetries);
    }
    return true;
  }
  else if((u[0]%based_nodes == u[1]%based_nodes && v[0]%based_nodes == v[1]%based_nodes) ||
          (u[0]%based_nodes == v[1]%based_nodes && v[0]%based_nodes == u[1]%based_nodes)){
    // A graph with symmetry can be created when using mutate_adjacency_1opt_s().
    // But I want mutate_adjacency_1opt_s() not to be called in this function
    // because I want the number of calls to 1opt_s() and 2opt_s() to be about the same.
    return false;
  }

  if(get_random(2) == 0){ // u[0]--v[1], v[0]--u[1]
    adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_VERTEX(v[1], u[0], nodes, symmetries);
    adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_VERTEX(v[0], u[1], nodes, symmetries);
    adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_VERTEX(u[1], v[0], nodes, symmetries);
    adjacency[v[1]%based_nodes][v_d[1]] = LOCAL_VERTEX(u[0], v[1], nodes, symmetries);
  }
  else{ // u[0]--u[1], v[0]--v[1]
    adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_VERTEX(u[1], u[0], nodes, symmetries);
    adjacency[u[1]%based_nodes][u_d[1]] = LOCAL_VERTEX(u[0], u[1], nodes, symmetries);
    adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_VERTEX(v[1], v[0], nodes, symmetries);
    adjacency[v[1]%based_nodes][v_d[1]] = LOCAL_VERTEX(v[0], v[1], nodes, symmetries);
  }

  if(check_isolated_vertex(u[0]%based_nodes, degree, adjacency) ||
     check_isolated_vertex(v[0]%based_nodes, degree, adjacency) ||
     check_isolated_vertex(u[1]%based_nodes, degree, adjacency) ||
     check_isolated_vertex(v[1]%based_nodes, degree, adjacency)){
    restore_adjacency(degree, (int *)adjacency);
    return false;
  }
  else{
    return true;
  }
}

static bool mutate_adjacency_1opt_s(const int nodes, const int degree,
                                    const int symmetries, int adjacency[nodes/symmetries][degree])
{
  int u[2], u_d[2], v[2], v_d[2]; // Declared with two elements since for backup_restore_adjacency()
  int based_nodes = nodes/symmetries;
  u[0]   = get_random(nodes);
  u_d[0] = get_random(degree);
  v[0]   = GLOBAL_ADJ(nodes, degree, symmetries, adjacency, u[0], u_d[0]);
  v_d[0] = get_degree_index(u[0], v[0], u_d[0], nodes, degree, symmetries, adjacency);

  backup_restore_adjacency(u, u_d, v, v_d, nodes, symmetries, MUTATE_1OPT);
  if(symmetries%2 == 0 && abs(u[0]-v[0]) == nodes/2) return false;

  int rnd   = (symmetries%2 == 1)? get_random(symmetries-1) : get_random(symmetries);
  int new_v = (rnd != symmetries-1)? v[0] + based_nodes*(rnd+1) : u[0] + based_nodes*(symmetries/2);
  //  int rnd   = get_random(symmetries-1);
  //  int new_v = v[0] + based_nodes*(rnd+1);
  int tmp_v = adjacency[u[0]%based_nodes][u_d[0]];
  int new_u = v[0] - (new_v - u[0]);
  int tmp_u = adjacency[v[0]%based_nodes][v_d[0]];

  if(v[0] == tmp_v && u[0] == tmp_u) return false; // No change
  
  adjacency[v[0]%based_nodes][v_d[0]] = LOCAL_VERTEX(new_u, v[0], nodes, symmetries);
  adjacency[u[0]%based_nodes][u_d[0]] = LOCAL_VERTEX(new_v, u[0], nodes, symmetries);
  
  return true;
}

void mutate_adjacency_general_s(const int nodes, const int degree,
				const int symmetries, int adjacency[nodes/symmetries][degree])
{
  CHECK_SYMMETRIES(nodes, symmetries);

  while(1){
    if(symmetries == 1){
      if(mutate_adjacency_2opt_s(nodes, degree, symmetries, adjacency))
	break;
    }
    else{
      if(get_random(2) == 0){
	if(mutate_adjacency_1opt_s(nodes, degree, symmetries, adjacency)){
	  break;
	}
      }
      else{
	if(mutate_adjacency_2opt_s(nodes, degree, symmetries, adjacency)){
	  break;
	}
      } 
    }
  }
}
