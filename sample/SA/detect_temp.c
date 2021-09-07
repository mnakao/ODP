#include "common.h"

double calc_min_temp_s()
{
  return -2.0 / log(0.0001);
}

double calc_min_temp()
{
  return calc_min_temp_s();
}

static bool _accept_s(const int nodes, const int current_diameter, const int diameter,
                      const double current_ASPL, const double ASPL, const int symmetries, double *max_diff_energy)
{
  if(diameter < current_diameter){
    return true;
  }
  else if(diameter > current_diameter){
    return false;
  }

  //  diameter == current_diameter
  double diff = (current_ASPL-ASPL)*nodes*(nodes-1)/symmetries;
  *max_diff_energy = MAX(*max_diff_energy, -1.0 * diff);

  return (ASPL <= current_ASPL);
}

static bool _accept(const int nodes, const int current_diameter, const int diameter,
                    const double current_ASPL, const double ASPL, double *max_diff_energy)
{
  return _accept_s(nodes, current_diameter, diameter, current_ASPL, ASPL, 1, max_diff_energy);
}

double calc_max_temp_s(const int nodes, const int degree, const int seed, const int symmetries)
{
  int lines = (nodes * degree)/2, diameter, current_diameter, ncalcs = 100;
  long sum;
  double ASPL, current_ASPL, max_diff_energy = 0;
  ODP_Restore r;

  ODP_Srand(seed);
  int (*edge)[2] = malloc(sizeof(int) * lines * 2);
  ODP_Generate_random_general_s(nodes, degree, symmetries, edge);
  int based_nodes = nodes/symmetries;
  int (*adjacency)[degree] = malloc(sizeof(int) * based_nodes * degree);

  ODP_Conv_edge2adjacency_general_s(nodes, lines, degree, edge, symmetries, adjacency);
  ODP_Init_aspl_general_s(nodes, degree, NULL, symmetries);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  current_diameter = diameter;
  current_ASPL     = ASPL;
  
  char *val = getenv("ODP_PROFILE");
  if(val) unsetenv("ODP_PROFILE");

  for(int i=0;i<ncalcs;i++){
    ODP_Mutate_adjacency_general_s(nodes, degree, NULL, symmetries, &r, adjacency);
    ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
    if(_accept_s(nodes, current_diameter, diameter, current_ASPL, ASPL, symmetries, &max_diff_energy)){
      current_diameter = diameter;
      current_ASPL     = ASPL;
    }
    else{
      ODP_Restore_adjacency_general(r, adjacency);
    }
  }
  ODP_Finalize_aspl();

  if(val) setenv("ODP_PROFILE", val, 1);

  free(edge);
  free(adjacency);
  return (-1.0 * max_diff_energy) / log(0.5);
}

double calc_max_temp(const int nodes, const int degree, const int seed)
{
  return calc_max_temp_s(nodes, degree, seed, 1);
}
