#include <mpi.h>
#include "common.h"
extern double calc_max_temp_mpi(int nodes, int degree, int seed, MPI_Comm comm);
extern double calc_min_temp_mpi();

static void print_help(char *argv)
{
  ERROR("%s [-N nodes] [-D degree] [-f input] [-o output] [-s seed] [-n calcs] [-w max_temp] [-c min_temp] [-Y] [-A]\n", argv);
}

static void set_args(const int argc, char **argv, int *nodes, int *degree, char **infname, char **outfname,
		     int *seed, long *ncalcs, double *max_temp, double *min_temp, bool *hill_climbing, bool *ASPL_priority)
{
  int result;
  while((result = getopt(argc,argv,"N:D:f:o:s:n:w:c:YA"))!=-1){
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
  int nodes = NOT_DEFINED, degree = NOT_DEFINED, lines, (*edge)[2], rank;
  int seed = 0, diameter, current_diameter, best_diameter, low_diameter;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = NOT_DEFINED, min_temp = NOT_DEFINED, ASPL, current_ASPL, best_ASPL, low_ASPL;
  ODP_Restore r;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  set_args(argc, argv, &nodes, &degree, &infname, &outfname, &seed,
	   &ncalcs, &max_temp, &min_temp, &hill_climbing, &ASPL_priority);
  
  ODP_Srand(seed);
  if(infname){
    if(nodes != NOT_DEFINED || degree != NOT_DEFINED)
      ERROR("When using -f option, you cannot use -N and -D.\n");
    lines = ODP_Get_lines(infname);
    edge = malloc(sizeof(int)*lines*2); // int edge[lines][2];
    ODP_Read_edge_general(infname, edge);
    nodes  = ODP_Get_nodes(lines, edge);
    degree = ODP_Get_degree(nodes, lines, edge);
  }
  else{
    if(nodes == NOT_DEFINED || degree == NOT_DEFINED)
      print_help(argv[0]);
    else if(nodes%2 == 1 && degree%2 == 1)
      ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);

    lines = (nodes * degree)/2;
    edge = malloc(sizeof(int)*lines*2); // int edge[lines][2];
    ODP_Generate_random_general(nodes, degree, edge);
  }
  
  if(max_temp == NOT_DEFINED)
    max_temp = calc_max_temp_mpi(nodes, degree, seed, MPI_COMM_WORLD);

  if(min_temp == NOT_DEFINED)
    min_temp = calc_min_temp_mpi();
  
  PRINT_R0("Nodes = %d, Degrees = %d\n", nodes, degree);
  PRINT_R0("Random seed = %d\n", seed);
  PRINT_R0("Number of calculations = %ld\n", ncalcs);
  if(hill_climbing)
    PRINT_R0("Method = Hill Climbing\n");
  else{
    PRINT_R0("Method = Simulated Annealing\n");
    PRINT_R0("Max, Min temperature = %f, %f\n", max_temp, min_temp);
  }
  if(infname)
    PRINT_R0("Input file name = %s\n", infname);
  if(outfname)
    PRINT_R0("Output file name = %s\n", outfname);
  
  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree);      // int adjacency[nodes][degree];
  int (*best_adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int best_adjacency[nodes][degree];

  ODP_Conv_edge2adjacency_general(nodes, lines, degree, edge, adjacency);
  ODP_Init_aspl_mpi_general(nodes, degree, NULL, MPI_COMM_WORLD);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);

  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL     = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);

  ODP_Set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
  MPI_Barrier(MPI_COMM_WORLD); // To measure time with ODP_PROFILE=1.
  double sa_time = get_time();
  if(diameter == low_diameter && ASPL == low_ASPL){
    PRINT_R0("Find optimum solution\n");
  }
  else{
    double cooling_rate = (hill_climbing)? 0 : pow(min_temp/max_temp, (double)1.0/ncalcs);
    double temp = max_temp;
    int	interval = (ncalcs < 100)? 1 : ncalcs/100;
    if(hill_climbing)
      PRINT_R0("Ncalcs : Best ASPL Gap ( Dia. )\n");
    else
      PRINT_R0("Ncalcs : Temp : current ASPL Gap ( Dia. ) : Best ASPL Gap ( Dia. )\n");

    for(long i=0;i<ncalcs;i++){
      if(i%interval == 0){
        if(hill_climbing)
          PRINT_R0("%ld\t%f ( %d )\n", i, 
                 best_ASPL-low_ASPL, best_diameter-low_diameter);
        else
          PRINT_R0("%ld\t%f\t%f ( %d )\t%f ( %d )\n", i, temp,
                 current_ASPL-low_ASPL, current_diameter-low_diameter,
                 best_ASPL-low_ASPL, best_diameter-low_diameter);
      }

      ODP_Mutate_adjacency_general(nodes, degree, NULL, &r, adjacency);
      ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
      if(diameter < best_diameter || (diameter == best_diameter && ASPL < best_ASPL)){
	best_diameter = diameter;
	best_sum      = sum;
	best_ASPL     = ASPL;
	memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);
	if(diameter == low_diameter && ASPL == low_ASPL){
          PRINT_R0("Find optimum solution\n");
	  break;
	}
      }
      
      if(accept(nodes, current_diameter, diameter, current_ASPL, ASPL, temp, hill_climbing, ASPL_priority)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
      }
      else{
	ODP_Restore_adjacency_general(r, adjacency);
      }
      temp *= cooling_rate;
    }
  }
  sa_time = get_time() - sa_time;
  ODP_Finalize_aspl();
  ODP_Conv_adjacency2edge_general(nodes, degree, NULL, best_adjacency, edge);

  PRINT_R0("---\n");
  PRINT_R0("Diameter        = %d\n", best_diameter);
  PRINT_R0("Diameter Gap    = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  PRINT_R0("ASPL            = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)nodes*(nodes-1)/2);
  PRINT_R0("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  PRINT_R0("Time            = %f sec.\n", sa_time);
  PRINT_R0("ASPL priority?  = %s\n", (ASPL_priority)? "Yes" : "No");
    //  PRINT_R0("Loop ?          = %s\n", (ODP_Check_loop(lines, edge))? "Yes" : "No");
    //  PRINT_R0("Multiple Edges? = %s\n", (ODP_Check_multiple_edges(lines, edge))? "Yes" : "No");
  
  if(outfname && rank == 0)
    ODP_Write_edge_general(lines, edge, outfname);
  
  free(edge);
  free(adjacency);
  free(best_adjacency);
  
  MPI_Finalize();
  return 0;
}
