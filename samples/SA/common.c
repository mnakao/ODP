#include "common.h"
static int _u[2], _u_d[2], _v[2], _v_d[2], _rnd;
static bool check_length(const int v, const int w, const int height, const int length)
{
  int w0 = WIDTH(v,height);
  int h0 = HEIGHT(v,height);
  int w1 = WIDTH(w,height);
  int h1 = HEIGHT(w,height);
  int distance = abs(w0 - w1) + abs(h0 - h1);

  return (distance <= length);
}

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

bool accept(const int nodes, const int current_diameter, const int diameter,
	    const double current_ASPL, const double ASPL, const double temp)
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
      double diff = (current_ASPL-ASPL)*nodes*(nodes-1);
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

static bool check_isolated_vertex(const int x, const int degree, const int (*adjacency)[degree])
{
  for(int i=1;i<degree;i++)
    if(adjacency[x][0] != adjacency[x][i])
      return false;

  return true;
}

static int get_degree_index(const int u, const int v, const int u_d,
                            const int degree, const int (*adjacency)[degree])
{
  if(u == v){ // loop
    for(int i=0;i<degree;i++)
      if(adjacency[v][i] == u && i != u_d)
  	return i;
  }
  else{
    for(int i=0;i<degree;i++)
      if(adjacency[v][i] == u)
    	return i;
  }

  ERROR("Something wrong !\n");
  return -1; // dummy
}

void restore_adjacency(const int degree, int *adjacency)
{
  adjacency[_u[0] * degree + _u_d[0]] = _v[0];
  adjacency[_v[0] * degree + _v_d[0]] = _u[0];
  adjacency[_u[1] * degree + _u_d[1]] = _v[1];
  adjacency[_v[1] * degree + _v_d[1]] = _u[1];
}

void mutate_adjacency_general(const int nodes, const int degree, int adjacency[nodes][degree])
{
  int u[2], v[2], u_d[2], v_d[2];
  while(1){
    while(1){
      u[0] = get_random(nodes);
      u[1] = get_random(nodes);
      if(u[0] == u[1]) continue;

      u_d[0] = get_random(degree);
      v[0] = adjacency[u[0]][u_d[0]];
      if(v[0] == u[1]) continue;

      u_d[1] = get_random(degree);
      v[1] = adjacency[u[1]][u_d[1]];
      if(v[1] == u[0] || v[0] == v[1]) continue;
      break;
    }

    v_d[0] = get_degree_index(u[0], v[0], u_d[0], degree, adjacency);
    v_d[1] = get_degree_index(u[1], v[1], u_d[1], degree, adjacency);

    // backup for restore_adjacency()
    memcpy(_u,   u, sizeof(int)*2);
    memcpy(_v,   v, sizeof(int)*2);
    memcpy(_u_d, u_d, sizeof(int)*2);
    memcpy(_v_d, v_d, sizeof(int)*2);
    _rnd = get_random(2);

    if(_rnd == 0){ // u[0]--v[1], v[0]--u[1]
      adjacency[u[0]][u_d[0]] = v[1];
      adjacency[u[1]][u_d[1]] = v[0];
      adjacency[v[0]][v_d[0]] = u[1];
      adjacency[v[1]][v_d[1]] = u[0];
    }
    else{ // u[0]--u[1], v[0]--v[1]
      adjacency[u[0]][u_d[0]] = u[1];
      adjacency[u[1]][u_d[1]] = u[0];
      adjacency[v[0]][v_d[0]] = v[1];
      adjacency[v[1]][v_d[1]] = v[0];
    }

    if(check_isolated_vertex(u[0], degree, adjacency) ||
       check_isolated_vertex(v[0], degree, adjacency) ||
       check_isolated_vertex(u[1], degree, adjacency) ||
       check_isolated_vertex(v[1], degree, adjacency)){
      restore_adjacency(degree, (int *)adjacency);
      continue;
    }
    break;
  }
}

void mutate_adjacency_grid(const int width, const int height, const int degree,
			   const int length, int (*adjacency)[degree])
{
  int nodes = width * height;
  while(1){
    mutate_adjacency_general(nodes, degree, adjacency);
    if(_rnd == 0 && check_length(_u[0], _v[1], height, length) && check_length(_u[1], _v[0], height, length)){
      break;
    }
    else if(_rnd == 1 && check_length(_u[0], _u[1], height, length) && check_length(_v[0], _v[1], height, length)){
      break;
    }
    else{
      restore_adjacency(degree, (int *)adjacency);
    }
  }
}
