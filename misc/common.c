#include "common.h"
#define MAX_NOPT 3
static int _u[MAX_NOPT], _v[MAX_NOPT], _u_d[MAX_NOPT], _v_d[MAX_NOPT];

void undo_adjacency(const int num, const int degree, int (*adjacency)[degree])
{
  for(int i=0;i<num;i++){
    adjacency[_u[i]][_u_d[i]] = _v[i];
    adjacency[_v[i]][_v_d[i]] = _u[i];
  }
}

void backup_adjacency(const int num, const int u[num], const int u_d[num], const int v[num], const int v_d[num])
{
  for(int i=0;i<num;i++){
    _u[i]   = u[i];
    _v[i]   = v[i];
    _u_d[i] = u_d[i];
    _v_d[i] = v_d[i];
  }
}

bool check_multiple_edges(const int u, const int u_d, const int v, 
                          const int degree, const int (*adjacency)[degree])
{
  for(int i=0;i<degree;i++)
    if(i!=u_d && adjacency[u][i] == v)
      return false;

  return true;
}

int get_degree_index(const int u, const int v, const int u_d, const int nodes,
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

  ERROR("Something Wrong ! [id=4]\n");
  return -1; // dummy
}

int get_random(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
}

double get_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1.0e-6 * t.tv_usec;
}

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

bool accept_s(const int nodes, const int current_diameter, const int diameter,
	      const double current_ASPL, const double ASPL, const double temp,
	      const bool hill_climbing, const bool ASPL_priority, const int symmetries)
{
  if(hill_climbing)
    return (ASPL <= current_ASPL);

  if(diameter < current_diameter && !ASPL_priority){
    return true;
  }
  else if(diameter > current_diameter && !ASPL_priority){
    return false;
  }
  else{ //  diameter == current_diameter
    if(ASPL <= current_ASPL){
      return true;
    }
    else{
      double diff = ((current_ASPL-ASPL)*nodes*(nodes-1))/symmetries;
      return exp(diff/temp) > uniform_rand();
    }
  }
}

bool accept(const int nodes, const int current_diameter, const int diameter, const double current_ASPL,
            const double ASPL, const double temp, const bool hill_climbing, const bool ASPL_priority)
{
  return accept_s(nodes, current_diameter, diameter, current_ASPL, ASPL, temp, hill_climbing, ASPL_priority, 1);
}  

bool accept_temp(const int nodes, const int current_diameter, const int diameter,
                 const double current_ASPL, const double ASPL, const double temp, double *max_diff_energy)
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
      double diff = ((current_ASPL-ASPL)*nodes*(nodes-1));
      *max_diff_energy = MAX(*max_diff_energy, -1.0 * diff);
       
      if(exp(diff/temp) > uniform_rand()){
        return true;
      }
      else{
        return false;
      }
    }
  }
}
