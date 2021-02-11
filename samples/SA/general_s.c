#include "common.h"

static void print_help(char *argv)
{
  ERROR("%s -N nodes -D degree -S symmetries [-o <output>] [-s <seed>] [-n <calcs>] [-w <max_temp>] [-c <min_temp>] [-A]\n", argv);
}

static void set_args(const int argc, char **argv, int *nodes, int *degree, int *symmetries, char *outfname, bool *enable_output,
		     int *seed, long *ncalcs, double *max_temp, double *min_temp, bool *enable_ASPL_priority)
{
  int result;
  while((result = getopt(argc,argv,"N:D:S:o:s:n:w:c:A"))!=-1){
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
      strcpy(outfname, optarg);
      *enable_output = true;
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
    case 'A':
      *enable_ASPL_priority = true;
      break;
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  char outfname[MAX_FILENAME_LENGTH];
  bool enable_ASPL_priority = false, enable_output = false;
  int nodes = NOT_DEFINED, degree = NOT_DEFINED, symmetries = NOT_DEFINED;
  int seed = 0, diameter, current_diameter, best_diameter, low_diameter;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = 100, min_temp = 0.22, ASPL, current_ASPL, best_ASPL, low_ASPL;

  set_args(argc, argv, &nodes, &degree, &symmetries, outfname, &enable_output,
	   &seed, &ncalcs, &max_temp, &min_temp, &enable_ASPL_priority);
  if(nodes == NOT_DEFINED || degree == NOT_DEFINED || symmetries == NOT_DEFINED)
    print_help(argv[0]);
  else if(nodes%2 == 1 && degree%2 == 1)
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
  ODP_Srand(seed);
  ODP_Generate_random_general_s(nodes, degree, symmetries, edge);
  create_time = get_time() - create_time;
  ODP_Conv_edge2adjacency_general_s(nodes, lines, degree, edge, symmetries, adjacency);

  ODP_Init_aspl_general_s(nodes, degree, NULL, symmetries); 
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL); 
  ODP_Conv_adjacency2edge_general_s(nodes, degree, NULL, adjacency, symmetries, edge);

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
    double cooling_rate = pow(min_temp/max_temp, (double)1.0/ncalcs);
    double temp = max_temp;
    int	interval = (ncalcs < 100)? 1 : ncalcs/100;
    printf("Ncalcs : Temp : Diameter Gap : ASPL Gap\n");
    for(long i=0;i<ncalcs;i++){
      if(i%interval == 0)
	printf("%ld\t%f\t%d\t%f\n", i, temp, best_diameter-low_diameter, best_ASPL-low_ASPL);

      ODP_Mutate_adjacency_general_s(nodes, degree, NULL, symmetries, adjacency);
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
      
      if(accept_s(nodes, current_diameter, diameter, current_ASPL, ASPL, temp, enable_ASPL_priority, symmetries)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
      }
      else{
	ODP_Restore_adjacency_general(adjacency);
      }
      temp *= cooling_rate;
    }
  }
  sa_time = get_time() - sa_time;
  ODP_Finalize_aspl();
  ODP_Conv_adjacency2edge_general_s(nodes, degree, NULL, best_adjacency, symmetries, edge);
  
  printf("---\n");
  printf("Diameter        = %d\n", best_diameter);
  printf("Diameter Gap    = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL            = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("Time            = %f/%f sec. (Create Graph/SA)\n", create_time, sa_time);
  printf("ASPL priority?  = %s\n", (enable_ASPL_priority)? "Yes" : "No");
  printf("Loop ?          = %s\n", (ODP_Check_loop(lines, edge))? "Yes" : "No");
  printf("Multiple Edges? = %s\n", (ODP_Check_multiple_edges(lines, edge))? "Yes" : "No");

  if(enable_output){
    ODP_Write_edge_general(lines, edge, outfname);
    printf("Generate ./%s\n", outfname);
  }
  
  free(edge);
  free(adjacency);
  free(best_adjacency);

  return 0;
}
