#include "common.h"

static void print_help(char *argv)
{
  ERROR("%s -N nodes -D degree [-o <output>] [-s <seed>] [-n <calcs>] [-w <max_temp>] [-c <min_temp>] [-A]\n", argv);
}

static void set_args(const int argc, char **argv, int *nodes, int *degree, char *fname,
		     int *seed, long *ncalcs, double *max_temp, double *min_temp, bool *enable_ASPL_priority)
{
  if(argc < 5)
    print_help(argv[0]);

  int result;
  while((result = getopt(argc,argv,"N:D:o:s:n:w:c:A"))!=-1){
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
  char *fname="general.edges";
  bool enable_ASPL_priority = false;
  int nodes, degree, seed = 0, diameter, current_diameter, best_diameter, low_diameter;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = 238.91, min_temp = 0.22, ASPL, current_ASPL, best_ASPL, low_ASPL;

  set_args(argc, argv, &nodes, &degree, fname, &seed, &ncalcs, &max_temp, &min_temp, &enable_ASPL_priority);
  if(nodes%2 == 1 && degree%2 == 1)
    ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);
  
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  printf("Max, Min temperature = %f, %f\n", max_temp, min_temp);
  
  int lines = (nodes * degree)/2;
  int (*edge)[2] = malloc(sizeof(int)*lines*2);                         // int edge[lines][2];
  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree);      // int adjacency[nodes][degree];
  int (*best_adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int best_adjacency[nodes][degree];

  double create_time = get_time();
  ODP_Generate_random_general(nodes, degree, seed, edge);
  create_time = get_time() - create_time;
  ODP_Conv_edge2adjacency(nodes, lines, degree, edge, adjacency);

  ODP_Init_aspl(nodes, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL     = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);

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

      ODP_mutate_adjacency_general(nodes, degree, NULL, adjacency);
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
      
      if(accept(nodes, current_diameter, diameter, current_ASPL, ASPL, temp, enable_ASPL_priority)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
      }
      else{
	ODP_restore_adjacency_general(nodes, degree, adjacency);
      }
      temp *= cooling_rate;
    }
  }
  sa_time = get_time() - sa_time;
  ODP_Finalize_aspl();
  ODP_Conv_adjacency2edge(nodes, degree, NULL, best_adjacency, edge);
  
  printf("---\n");
  printf("Diameter       = %d\n", best_diameter);
  printf("Diameter Gap   = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL           = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap       = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("Time           = %f/%f sec. (Create Graph/SA)\n", create_time, sa_time);
  printf("ASPL priority? = %s\n", (enable_ASPL_priority)? "Yes" : "No");

  //  ODP_Write_edge_general(lines, edge, fname);
  //  printf("Generate ./%s\n", fname);
  
  free(edge);
  free(adjacency);
  free(best_adjacency);

  return 0;
}
