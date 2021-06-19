#include "common.h"

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
