#include "common.h"

static void print_help(char *argv)
{
  ERROR("%s [-N nodes] [-D degree] [-s seed] [-n calcs]\n", argv);
}

static void set_args(const int argc, char **argv, int *nodes, int *degree,int *seed, long *ncalcs)
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
  int seed = 0, diameter, low_diameter;
  long sum, ncalcs = 10000;
  double ASPL, new_ASPL, low_ASPL, max_diff_energy = 0;

  set_args(argc, argv, &nodes, &degree, &seed, &ncalcs);
  
  if(nodes == NOT_DEFINED || degree == NOT_DEFINED)
    print_help(argv[0]);
  else if(nodes%2 == 1 && degree%2 == 1)
    ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);
  
  lines = (nodes * degree)/2;
  edge  = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ODP_Srand(seed);
  ODP_Generate_random_general(nodes, degree, edge);
  
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);

  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree);      // int adjacency[nodes][degree];
  ODP_Conv_edge2adjacency_general(nodes, lines, degree, edge, adjacency);
  ODP_Init_aspl_general(nodes, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  ODP_Set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
  
  if(diameter == low_diameter && ASPL == low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    for(long i=0;i<ncalcs;i++){
      ODP_Mutate_adjacency_general(nodes, degree, NULL, NULL, adjacency);
      ODP_Set_aspl(adjacency, &diameter, &sum, &new_ASPL);
      double diff = (new_ASPL - ASPL) * nodes * (nodes-1);
      max_diff_energy = MAX(max_diff_energy, diff);
      new_ASPL = ASPL;
    }
  }
  ODP_Finalize_aspl();
  
  printf("---\n");
  // Set max temperature to accept it   50% in maximum diff energy.
  printf("Proposed max temperature is %f\n", (-1.0 * max_diff_energy) / log(0.5));
  // Set min temperature to accept it 0.01% in minimum diff energy.
  printf("Proposed min temperature is %f\n", (-2.0) / log(0.0001));
  
  free(edge);
  free(adjacency);

  return 0;
}
