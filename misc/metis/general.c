#include "common.h"
#include "metis.h"
extern double calc_max_temp(int nodes, int degree, int seed);
extern double calc_min_temp();

static void print_help(char *argv)
{
  ERROR("%s [-N nodes] [-D degree] [-f input] [-o output] [-s seed] [-n calcs] [-w max_temp] [-c min_temp]\n", argv);
}

static void set_args(const int argc, char **argv, int *nodes, int *degree, char **infname, char **outfname,
		     int *seed, long *ncalcs, double *max_temp, double *min_temp)
{
  int result;
  while((result = getopt(argc,argv,"N:D:f:o:s:n:w:c:"))!=-1){
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
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  char *infname = NULL, *outfname = NULL;
  int nodes = NOT_DEFINED, degree = NOT_DEFINED, lines, (*edge)[2];
  int seed = 0, diameter, current_diameter, best_diameter, low_diameter, ncuts, best_ncuts, current_ncuts;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = NOT_DEFINED, min_temp = NOT_DEFINED, ASPL, current_ASPL, best_ASPL, low_ASPL;
  ODP_Restore r;

  set_args(argc, argv, &nodes, &degree, &infname, &outfname,
           &seed, &ncalcs, &max_temp, &min_temp);
  
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

  if(sizeof(int) != sizeof(idx_t)) exit(1);
  if(max_temp == NOT_DEFINED)
    max_temp = calc_max_temp(nodes, degree, seed);
  
  if(min_temp == NOT_DEFINED)
    min_temp = calc_min_temp();
  
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  printf("Max, Min temperature = %f, %f\n", max_temp, min_temp);

  if(infname)
    printf("Input file name = %s\n", infname);
  if(outfname)
    printf("Output file name = %s\n", outfname);

  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree);      // int adjacency[nodes][degree];
  int (*best_adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int best_adjacency[nodes][degree];

  ODP_Conv_edge2adjacency_general(nodes, lines, degree, edge, adjacency);
  ODP_Init_aspl_general(nodes, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);

  idx_t nvtxs = nodes, ncon = 1, nparts = 2, part[nodes];
  idx_t *xadj = malloc(sizeof(idx_t) * (nodes + 1));
  for(int i=0;i<nodes+1;i++)
    xadj[i] = i * degree;

  METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, (idx_t*)adjacency, NULL, NULL, NULL,
                           &nparts, NULL, NULL, NULL, &ncuts, part);

  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL     = ASPL;
  best_ncuts    = current_ncuts    = ncuts;
  memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);

  ODP_Set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);

  double sa_time = get_time();
  if(diameter == low_diameter && ASPL == low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    double cooling_rate = pow(min_temp/max_temp, (double)1.0/ncalcs);
    double temp = max_temp;
    int	interval = (ncalcs < 100)? 1 : ncalcs/100;
    printf("Ncalcs : Temp : current ASPL Gap ( Dia. ) : Best ASPL Gap ( Dia. ) : current ncuts : best ncuts\n");
    for(long i=0;i<ncalcs;i++){
      if(i%interval == 0){
        printf("%ld\t%f\t%f ( %d )\t%f ( %d ) %d %d\n", i, temp,
               current_ASPL-low_ASPL, current_diameter-low_diameter,
               best_ASPL-low_ASPL, best_diameter-low_diameter, current_ncuts, best_ncuts);
      }

      ODP_Srand(seed+i);
      ODP_Mutate_adjacency_general(nodes, degree, NULL, &r, adjacency);
      ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
      METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, (idx_t*)adjacency, NULL, NULL, NULL,
                               &nparts, NULL, NULL, NULL, &ncuts, part);

      if(diameter < best_diameter || (diameter == best_diameter && ncuts > best_ncuts)){
	best_diameter = diameter;
	best_sum      = sum;
	best_ASPL     = ASPL;
        best_ncuts    = ncuts;
	memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);
	if(diameter == low_diameter && ASPL == low_ASPL){
	  printf("Find optimum solution\n");
	  break;
	}
      }
      
      if(accept(nodes, current_diameter, diameter, current_ncuts, ncuts, temp)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
        current_ncuts    = ncuts;
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
  
  printf("---\n");
  printf("Diameter        = %d\n", best_diameter);
  printf("Diameter Gap    = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL            = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("ncuts           = %d\n", best_ncuts);
  printf("Time            = %f sec.\n", sa_time);

  if(outfname)
    ODP_Write_edge_general(lines, edge, outfname);
  
  free(edge);
  free(adjacency);
  free(best_adjacency);

  return 0;
}
