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
	      const int current_ncuts, const int ncuts, const double temp, const int symmetries)
{
  if(diameter < current_diameter){
    return true;
  }
  else if(diameter > current_diameter){
    return false;
  }
  else{ //  diameter == current_diameter
    if(ncuts >= current_ncuts){
      return true;
    }
    else{
      double diff = ((double)ncuts-current_ncuts)*sqrt(nodes);
      return exp(diff/temp) > uniform_rand();
    }
  }
}

bool accept(const int nodes, const int current_diameter, const int diameter, const int current_ncuts,
            const int ncuts, const double temp)
{
  return accept_s(nodes, current_diameter, diameter, current_ncuts, ncuts, temp, 1);
}  
