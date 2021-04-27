#include "common.h"

static void print_help(char *argv)
{
  ERROR("%s [-N nodes] [-D degree] [-s seed] [-n calcs]\n", argv);
}

static void set_args(const int argc, char **argv, int *nodes, int *degree, int *seed, long *ncalcs)
{
  int result;
  while((result = getopt(argc,argv,"N:D:s:n:"))!=-1){
    switch(result){
    case 'N':
      *nodes = atoi(optarg);
      if(*nodes <= 0)
        ERROR("-N value > 0\n");
      break;
    case 'D':
      *degree = atoi(optarg);
      if(*degree <= 0)
        ERROR("-D value > 0\n");
      break;
    case 's':
      *seed = atoi(optarg);
      if(*seed < 0)
        ERROR("-s value >= 0\n");
      break;
    case 'n':
      *ncalcs = atol(optarg);
      if(*ncalcs < 0)
        ERROR("-n value >= 0\n");
      break;
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  int nodes = NOT_DEFINED, degree = NOT_DEFINED, lines, (*edge)[2];
  int seed = 0, diameter, current_diameter, best_diameter;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = 100, min_temp = 0.22, ASPL, current_ASPL, best_ASPL, max_diff_energy = 0;
  ODP_Restore r;

  set_args(argc, argv, &nodes, &degree, &seed, &ncalcs);
  
  ODP_Srand(seed);
  if(nodes == NOT_DEFINED || degree == NOT_DEFINED)
    print_help(argv[0]);
  else if(nodes%2 == 1 && degree%2 == 1)
    ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);

  lines = (nodes * degree)/2;
  edge = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ODP_Generate_random_general(nodes, degree, edge);
  
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);

  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree);      // int adjacency[nodes][degree];
  int (*best_adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int best_adjacency[nodes][degree];

  ODP_Conv_edge2adjacency_general(nodes, lines, degree, edge, adjacency);
  ODP_Init_aspl_general(nodes, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);

  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL     = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);

  double cooling_rate = pow(min_temp/max_temp, (double)1.0/ncalcs);
  double temp = max_temp;
  for(long i=0;i<ncalcs;i++){
    ODP_Mutate_adjacency_general(nodes, degree, NULL, &r, adjacency);
    ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
    if(diameter < best_diameter || (diameter == best_diameter && ASPL < best_ASPL)){
      best_diameter = diameter;
      best_sum      = sum;
      best_ASPL     = ASPL;
      memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);
    }
    
    if(accept_temp(nodes, current_diameter, diameter, current_ASPL, ASPL, temp, &max_diff_energy)){
      current_diameter = diameter;
      current_ASPL     = ASPL;
    }
    else{
      ODP_Restore_adjacency_general(r, adjacency);
    }
    temp *= cooling_rate;
  }
  ODP_Finalize_aspl();

  // Set max temperature to accept it   50% in maximum diff energy.
  printf("Proposed max temperature is %f\n", (-1.0 * max_diff_energy) / log(0.5));
  // Set min temperature to accept it 0.01% in minimum diff energy.
  printf("Proposed min temperature is %f\n", (-2.0) / log(0.0001));
    
  free(edge);
  free(adjacency);
  free(best_adjacency);

  return 0;
}
