#include "common_s.h"

static void print_help(char *argv)
{
  ERROR("%s -N nodes -D degree -S symmetries [-o <output>] [-s <seed>] [-n <calcs>] [-w <max_temp>] [-c <min_temp>]\n", argv);
}

static void set_args(const int argc, char **argv, int *nodes, int *degree, int *symmetries,
		     char *fname, int *seed, long *ncalcs, double *max_temp, double *min_temp)
{
  if(argc < 7)
    print_help(argv[0]);

  int result;
  while((result = getopt(argc,argv,"N:D:S:o:s:n:w:c:"))!=-1){
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
    case 'S':
      *symmetries = atoi(optarg);
      if(*symmetries <= 0)
        ERROR("-S value > 0\n");
      break;
    case 'o':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Output filename is long (%s).\n", optarg);
      strcpy(fname, optarg);
      break;
    case 's':
      *seed = atoi(optarg);
      if(*seed < 0)
        ERROR("-s value >= 0\n");
      break;
    case 'n':
      *ncalcs = atol(optarg);
      if(*ncalcs <= 0)
        ERROR("-n value >= 1\n");
      break;
    case 'w':
      *max_temp = atof(optarg);
      if(*max_temp <= 0)
        ERROR("-w value > 0\n");
      break;
    case 'c':
      *min_temp = atof(optarg);
      if(*min_temp <= 0)
        ERROR("MIN value > 0\n");
      break;
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  char *fname="general_s.edges";
  int nodes, degree, symmetries, seed = 0, diameter, current_diameter, best_diameter, low_diameter;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = 238.91, min_temp = 0.217147, ASPL, current_ASPL, best_ASPL, low_ASPL;

  set_args(argc, argv, &nodes, &degree, &symmetries, fname, &seed, &ncalcs, &max_temp, &min_temp);
  if(nodes%2 == 1 && degree%2 == 1)
    ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);
  else if(nodes%symmetries != 0)
    ERROR("Invalid nodes(%d) or symmetries(%d)\n", nodes, symmetries);

  printf("Nodes = %d, Degrees = %d, Symmetries = %d\n", nodes, degree, symmetries);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  printf("Max, Min temperature = %f, %f\n", max_temp, min_temp);
  
  int lines = (nodes * degree)/2;
  int based_nodes = nodes/symmetries;
  int (*edge)[2] = malloc(sizeof(int)*lines*2);                           // int edge[lines][2];
  int (*adjacency)[degree] = malloc(sizeof(int)*based_nodes*degree);      // int adjacency[based_nodes][degree];
  int (*best_adjacency)[degree] = malloc(sizeof(int)*based_nodes*degree); // int best_adjacency[based_nodes][degree];

  double create_time = get_time();
  ODP_Generate_random_general_s(nodes, degree, seed, symmetries, edge);
  create_time = get_time() - create_time;
  ODP_Conv_edge2adjacency_s(nodes, lines, edge, symmetries, adjacency);

  ODP_Init_aspl_s(nodes, degree, NULL, symmetries);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  ODP_Conv_adjacency2edge_s(nodes, degree, NULL, adjacency, symmetries, edge);

  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL     = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int)*based_nodes*degree);

  ODP_Set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
  double sa_time = get_time();
  if(diameter == low_diameter && ASPL == low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    double cooling_rate = (max_temp != min_temp)? pow(min_temp/max_temp, (double)1.0/ncalcs) : 1.0;
    double temp = max_temp;
    printf("Ncalcs : Temp : Diameter Gap : ASPL Gap\n");
    for(long i=0;i<ncalcs;i++){
      if(i%10000 == 0)
	printf("%ld\t%f\t%d\t%f\n", i, temp, best_diameter-low_diameter, best_ASPL-low_ASPL);

      mutate_adjacency_general_s(nodes, degree, symmetries, adjacency);
      ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);

      if(diameter < best_diameter || (diameter == best_diameter && ASPL < best_ASPL)){
	best_diameter = diameter;
	best_sum      = sum;
	best_ASPL     = ASPL;
	memcpy(best_adjacency, adjacency, sizeof(int)*based_nodes*degree);
	if(diameter == low_diameter && ASPL == low_ASPL){
	  printf("Find optimum solution\n");
	  break;
	}
      }
      
      if(accept_s(nodes, current_diameter, diameter, current_ASPL, ASPL, temp, symmetries)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
      }
      else{
	restore_adjacency(degree, (int *)adjacency);
      }
      temp *= cooling_rate;
    }
  }
  sa_time = get_time() - sa_time;
  ODP_Finalize_aspl();
  ODP_Conv_adjacency2edge_s(nodes, degree, NULL, best_adjacency, symmetries, edge);
  
  printf("---\n");
  printf("Diameter       = %d\n", best_diameter);
  printf("Diameter Gap   = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL           = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap       = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("Time           = %f/%f sec. (Create Graph/SA)\n", create_time, sa_time);
  
  //  ODP_Write_edge_general(lines, edge, fname);
  //  printf("Generate ./%s\n", fname);
  
  free(edge);
  free(adjacency);
  free(best_adjacency);

  return 0;
}
