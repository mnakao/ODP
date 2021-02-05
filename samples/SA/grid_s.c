#include "common.h"

static void print_help(char *argv)
{
  ERROR("%s -W width -H height -D degree -L length -S symmetries [-o <output>] [-s <seed>] [-n <calcs>] [-w <max_temp>] [-c <min_temp>] [-A]\n", argv);
}

static void set_args(const int argc, char **argv, int *width, int *height, int *degree, int *length, int *symmetries,
		     char *fname, int *seed, long *ncalcs, double *max_temp, double *min_temp, bool *enable_ASPL_priority)
{
  int result;
  while((result = getopt(argc,argv,"W:H:D:L:S:o:s:n:w:c:A"))!=-1){
    switch(result){
    case 'W':
      *width = atoi(optarg);
      if(*width <= 0)
        ERROR("-W value > 0\n");
      break;
    case 'H':
      *height = atoi(optarg);
      if(*height <= 0)
        ERROR("-H value > 0\n");
      break;
    case 'D':
      *degree = atoi(optarg);
      if(*degree <= 0)
        ERROR("-D value > 0\n");
      break;
    case 'L':
      *length = atoi(optarg);
      if(*length <= 0)
        ERROR("-L value > 0\n");
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
  char *fname="grid_s.edges";
  bool enable_ASPL_priority = false;
  int width = NOT_DEFINED, height = NOT_DEFINED, degree = NOT_DEFINED, length = NOT_DEFINED, symmetries = NOT_DEFINED;
  int seed = 0, diameter, current_diameter, best_diameter, low_diameter;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = 100, min_temp = 0.22, ASPL, current_ASPL, best_ASPL, low_ASPL;

  set_args(argc, argv, &width, &height, &degree, &length, &symmetries, fname, &seed, &ncalcs, &max_temp, &min_temp, &enable_ASPL_priority);
  int nodes = width * height;
  if(width == NOT_DEFINED || height == NOT_DEFINED || degree == NOT_DEFINED || length == NOT_DEFINED || symmetries == NOT_DEFINED)
    print_help(argv[0]);
  else if(nodes%2 == 1 && degree%2 == 1)
    ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);
  else if(nodes%symmetries != 0)
    ERROR("Invalid nodes(%d) or symmetries(%d)\n", nodes, symmetries);
  
  printf("Width = %d, Height = %d, Degrees = %d, Length = %d, Symmetries = %d\n",
	 width, height, degree, length, symmetries);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  printf("Max, Min temperature = %f, %f\n", max_temp, min_temp);
  
  int lines = (nodes * degree)/2;
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  int based_nodes = nodes/symmetries;
  int (*adjacency)[degree] = malloc(sizeof(int) * based_nodes * degree); // int adjacency[based_nodes][degree];
  int (*best_adjacency)[degree] = malloc(sizeof(int) * based_nodes * degree); // int best_adjacency[based_nodes][degree];

  double create_time = get_time();
  ODP_Srand(seed);
  ODP_Generate_random_grid_s(width, height, degree, length, symmetries, edge);
  create_time = get_time() - create_time;
  ODP_Conv_edge2adjacency_grid_s(width, height, lines, degree, edge, symmetries, adjacency);
  
  ODP_Init_aspl(nodes, degree, NULL);
  
  int adjacency2[nodes][degree], edge2[lines][2];
  ODP_Conv_adjacency2edge_grid_s(width, height, degree, NULL, adjacency, symmetries, edge2);
  ODP_Conv_edge2adjacency(nodes, lines, degree, edge2, adjacency2);
  ODP_Set_aspl(adjacency2, &diameter, &sum, &ASPL);
  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL     = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int)*based_nodes*degree);

  ODP_Set_lbounds_grid(width, height, degree, length, &low_diameter, &low_ASPL);
  double sa_time = get_time();
  if(diameter == low_diameter && ASPL == low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    double cooling_rate = (max_temp != min_temp)? pow(min_temp/max_temp, (double)1.0/ncalcs) : 1.0;
    double temp = max_temp;
    printf("Ncalcs : Temp : Diameter Gap : ASPL Gap\n");
    for(long i=0;i<ncalcs;i++){
      if(i%(ncalcs/100) == 0)
	printf("%ld\t%f\t%d\t%f\n", i, temp, best_diameter-low_diameter, best_ASPL-low_ASPL);

      ODP_Mutate_adjacency_grid_s(width, height, degree, NULL, length, symmetries, adjacency);
      
      ODP_Conv_adjacency2edge_grid_s(width, height, degree, NULL, adjacency, symmetries, edge2);
      ODP_Conv_edge2adjacency(nodes, lines, degree, edge2, adjacency2);
      ODP_Set_aspl(adjacency2, &diameter, &sum, &ASPL);

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
	ODP_Restore_adjacency_grid_s(width, height, degree, symmetries, adjacency);
      }
      temp *= cooling_rate;
    }
  }
  sa_time = get_time() - sa_time;  
  ODP_Finalize_aspl();
  ODP_Conv_adjacency2edge_grid_s(width, height, degree, NULL, adjacency, symmetries, edge);
  
  printf("---\n");
  printf("Diameter        = %d\n", best_diameter);
  printf("Diameter Gap    = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL            = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("Time            = %f/%f sec. (Create Graph/SA)\n", create_time, sa_time);
  printf("ASPL Priority?  = %s\n", (enable_ASPL_priority)? "Yes" : "No");
  printf("Loop ?          = %s\n", (ODP_Check_loop(lines, edge))? "Yes" : "No");
  printf("Multiple Edges? = %s\n", (ODP_Check_multiple_edges(lines, edge))? "Yes" : "No");
  
  //  ODP_Write_edge_grid(lines, height, edge, fname);
  //  printf("Generate ./%s\n", fname);
  
  free(edge);
  free(adjacency);
  free(best_adjacency);

  return 0;
}
