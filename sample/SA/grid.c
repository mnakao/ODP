#include "common.h"
extern double calc_grid_max_temp(const int width, const int height, const int degree, const int length, const int seed);
extern double calc_grid_min_temp();

static void print_help(char *argv)
{
  ERROR("%s -L length [-W width] [-H height] [-D degree] [-f input] [-o output] [-s seed] [-n calcs] [-w max_temp] [-c min_temp] [-Y] [-A]\n", argv);
}

static void set_args(const int argc, char **argv, int *width, int *height, int *degree, int *length, char **infname, char **outfname,
		     int *seed, long *ncalcs, double *max_temp, double *min_temp, bool *hill_climbing, bool *ASPL_priority)
{
  int result;
  while((result = getopt(argc,argv,"W:H:D:L:f:o:s:n:w:c:YA"))!=-1){
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
    case 'Y':
      *hill_climbing = true;
      break;
    case 'A':
      *ASPL_priority = true;
      break;
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  char *infname = NULL, *outfname = NULL;
  bool hill_climbing = false, ASPL_priority = false;
  int width = NOT_DEFINED, height = NOT_DEFINED, degree = NOT_DEFINED, length = NOT_DEFINED, lines, (*edge)[2], nodes;
  int seed = 0, diameter, current_diameter, best_diameter, low_diameter;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = NOT_DEFINED, min_temp = NOT_DEFINED, ASPL, current_ASPL, best_ASPL, low_ASPL;
  ODP_Restore r;

  set_args(argc, argv, &width, &height, &degree, &length, &infname, &outfname, &seed,
           &ncalcs, &max_temp, &min_temp, &hill_climbing, &ASPL_priority);

  ODP_Srand(seed);
  if(infname){
    lines = ODP_Get_lines(infname);
    edge = malloc(sizeof(int)*lines*2); // int edge[lines][2];
    ODP_Read_edge_grid(infname, &width, &height, edge);
    nodes  = ODP_Get_nodes(lines, edge);
    degree = ODP_Get_degree(nodes, lines, edge);
  }
  else{
    nodes = width * height;
    if(width == NOT_DEFINED || height == NOT_DEFINED || degree == NOT_DEFINED || length == NOT_DEFINED)
      print_help(argv[0]);
    else if(nodes%2 == 1 && degree%2 == 1)
      ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);

    lines = (nodes * degree)/2;
    edge = malloc(sizeof(int)*lines*2); // int edge[lines][2];
    ODP_Generate_random_grid(width, height, degree, length, edge);
  }

  if(max_temp == NOT_DEFINED)
    max_temp = calc_grid_max_temp(width, height, degree, length, seed);
  
  if(min_temp == NOT_DEFINED)
    min_temp = calc_grid_min_temp();
  
  printf("Width = %d, Height = %d, Degrees = %d, Length = %d\n",
	 width, height, degree, length);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  if(hill_climbing)
    printf("Method = Hill Climbing\n");
  else{
    printf("Method = Simulated Annealing\n");
    printf("Max, Min temperature = %f, %f\n", max_temp, min_temp);
  }
  if(infname)
    printf("Input file name = %s\n", infname);
  if(outfname)
    printf("Output file name = %s\n", outfname);
  
  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int adjacency[nodes][degree];
  int (*best_adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int best_adjacency[nodes][degree];

  ODP_Conv_edge2adjacency_grid(width, height, lines, degree, edge, adjacency);
  ODP_Init_aspl_grid(width, height, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  
  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL     = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);

  ODP_Set_lbounds_grid(width, height, degree, length, &low_diameter, &low_ASPL);
  double sa_time = get_time();
  if(diameter == low_diameter && ASPL == low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    double cooling_rate = (hill_climbing)? 0 : pow(min_temp/max_temp, (double)1.0/ncalcs);
    double temp = max_temp;
    int interval = (ncalcs < 100)? 1 : ncalcs/100;
    if(hill_climbing)
      printf("Ncalcs : Best ASPL Gap ( Dia. )\n");
    else
      printf("Ncalcs : Temp : current ASPL Gap ( Dia. ) : Best ASPL Gap ( Dia. )\n");
    for(long i=0;i<ncalcs;i++){
      if(i%interval == 0){
        if(hill_climbing)
          printf("%ld\t%f ( %d )\n", i,
                 best_ASPL-low_ASPL, best_diameter-low_diameter);
        else
          printf("%ld\t%f\t%f ( %d )\t%f ( %d )\n", i, temp,
                 current_ASPL-low_ASPL, current_diameter-low_diameter,
                 best_ASPL-low_ASPL, best_diameter-low_diameter);
      }

      ODP_Mutate_adjacency_grid(width, height, degree, NULL, length, &r, adjacency);
      ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
      if(diameter < best_diameter || (diameter == best_diameter && ASPL < best_ASPL)){
	best_diameter = diameter;
	best_sum      = sum;
	best_ASPL     = ASPL;
	memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);
	if(diameter == low_diameter && ASPL == low_ASPL){
	  printf("Find optimum solution\n");
	  break;
	}
      }
      
      if(accept(nodes, current_diameter, diameter, current_ASPL, ASPL, temp, hill_climbing, ASPL_priority)){
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
  ODP_Conv_adjacency2edge_grid(width, height, degree, NULL, best_adjacency, edge);
  
  printf("---\n");
  printf("Diameter        = %d\n", best_diameter);
  printf("Diameter Gap    = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL            = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("Time            = %f sec.\n", sa_time);
  printf("ASPL Priority?  = %s\n", (ASPL_priority)? "Yes" : "No");
  //  printf("Loop ?          = %s\n", (ODP_Check_loop(lines, edge))? "Yes" : "No");
  //  printf("Multiple Edges? = %s\n", (ODP_Check_multiple_edges(lines, edge))? "Yes" : "No");

  if(outfname)
    ODP_Write_edge_grid(lines, height, edge, outfname);
  
  free(edge);
  free(adjacency);
  free(best_adjacency);

  return 0;
}
