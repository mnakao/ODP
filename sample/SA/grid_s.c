#include "common.h"

static void print_help(char *argv)
{
  ERROR("%s -L length -S symmetries [-W width] [-H height] [-D degree] [-o output] [-s seed] [-n calcs] [-w max_temp] [-c min_temp] [-A]\n", argv);
}

static void set_args(const int argc, char **argv, int *width, int *height, int *degree, int *length, int *symmetries, char **infname,
                     char **outfname, int *seed, long *ncalcs, double *max_temp, double *min_temp, bool *enable_ASPL_priority)
{
  int result;
  while((result = getopt(argc,argv,"W:H:D:L:S:f:o:s:n:w:c:A"))!=-1){
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
    case 'f':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Input filename is long (%s).\n", optarg);
      *infname = malloc(MAX_FILENAME_LENGTH);
      strcpy(*infname, optarg);
      break;
    case 'o':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Output filename is long (%s).\n", optarg);
      *outfname = malloc(MAX_FILENAME_LENGTH);
      strcpy(*outfname, optarg);
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
  char *infname = NULL, *outfname = NULL;
  bool enable_ASPL_priority = false;
  int width = NOT_DEFINED, height = NOT_DEFINED, degree = NOT_DEFINED, length = NOT_DEFINED, symmetries = 1, lines, (*edge)[2], nodes;
  int seed = 0, diameter, current_diameter, best_diameter, low_diameter;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = 100, min_temp = 0.22, ASPL, current_ASPL, best_ASPL, low_ASPL;
  ODP_Restore r;

  set_args(argc, argv, &width, &height, &degree, &length, &symmetries, &infname,
           &outfname, &seed, &ncalcs, &max_temp, &min_temp, &enable_ASPL_priority);
  
  ODP_Srand(seed);
  if(infname){
    lines = ODP_Get_lines(infname);
    edge = malloc(sizeof(int)*lines*2); // int edge[lines][2];
    ODP_Read_edge_grid(infname, &width, &height, edge);
    nodes  = ODP_Get_nodes(lines, edge);
    degree = ODP_Get_degree(nodes, lines, edge);
    if(nodes%symmetries != 0)
      ERROR("Invalid nodes(%d) or symmetries(%d)\n", nodes, symmetries);
  }
  else{
    nodes = width * height;
    if(width == NOT_DEFINED || height == NOT_DEFINED || degree == NOT_DEFINED || length == NOT_DEFINED)
      print_help(argv[0]);
    else if(nodes%2 == 1 && degree%2 == 1)
      ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);
    else if(nodes%symmetries != 0)
      ERROR("Invalid nodes(%d) or symmetries(%d)\n", nodes, symmetries);
    
    lines = (nodes * degree)/2;
    edge = malloc(sizeof(int)*lines*2); // int edge[lines][2];
    ODP_Generate_random_grid_s(width, height, degree, length, symmetries, edge);
  }
  
  printf("Width = %d, Height = %d, Degrees = %d, Length = %d, Symmetries = %d\n",
	 width, height, degree, length, symmetries);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  printf("Max, Min temperature = %f, %f\n", max_temp, min_temp);
  if(infname)
    printf("Input file name = %s\n", infname);
  if(outfname)
    printf("Output file name = %s\n", outfname);
  
  int based_nodes = nodes/symmetries;
  int (*adjacency)[degree] = malloc(sizeof(int) * based_nodes * degree); // int adjacency[based_nodes][degree];
  int (*best_adjacency)[degree] = malloc(sizeof(int) * based_nodes * degree); // int best_adjacency[based_nodes][degree];

  ODP_Conv_edge2adjacency_grid_s(width, height, lines, degree, edge, symmetries, adjacency);
  ODP_Init_aspl_grid_s(width, height, degree, NULL, symmetries);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);

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
    double cooling_rate = pow(min_temp/max_temp, (double)1.0/ncalcs);
    double temp = max_temp;
    int interval = (ncalcs < 100)? 1 : ncalcs/100;
    printf("Ncalcs : Temp : Diameter Gap : ASPL Gap\n");
    for(long i=0;i<ncalcs;i++){
      if(i%interval == 0)
	printf("%ld\t%f\t%d\t%f\n", i, temp, best_diameter-low_diameter, best_ASPL-low_ASPL);

      ODP_Mutate_adjacency_grid_s(width, height, degree, NULL, length, symmetries, &r, adjacency);
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
	ODP_Restore_adjacency_grid(r, adjacency);
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
  printf("Time            = %f sec.\n", sa_time);
  printf("ASPL Priority?  = %s\n", (enable_ASPL_priority)? "Yes" : "No");
  //  printf("Loop ?          = %s\n", (ODP_Check_loop(lines, edge))? "Yes" : "No");
  //  printf("Multiple Edges? = %s\n", (ODP_Check_multiple_edges(lines, edge))? "Yes" : "No");
  
  if(outfname)
    ODP_Write_edge_grid(lines, height, edge, outfname);
  
  free(edge);
  free(adjacency);
  free(best_adjacency);

  return 0;
}
