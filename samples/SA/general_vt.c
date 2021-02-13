#include "common.h"
#define USED     true
#define NOT_USED false
int _v, _d, _r;

static int get_random(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
}

static void restore_adjacency_general(const int nodes, const int degree, bool *used, int *adjacency)
{
  adjacency[_d] = _v;
  adjacency[_d+degree/2] = nodes - _v;
  used[_v]      = USED;
  used[_r]      = NOT_USED;
}

static void mutate_adjacency_general(const int nodes, const int degree, bool used[nodes/2], int adjacency[nodes])
{
  while(1){
    int d = get_random(degree/2);
    int r = get_random((nodes-1)/2) + 1;
    if(used[r] == NOT_USED){
      _v = adjacency[d];
      _d = d;
      _r = r;
      used[adjacency[d]]    = NOT_USED;
      used[r]               = USED;
      adjacency[d]          = r;
      adjacency[d+degree/2] = nodes - r;
      break;
    }
  }
}

static void init_adjacency_general(const int nodes, const int degree, bool used[nodes/2], int adjacency[nodes])
{
  if(degree%2 == 1)
    adjacency[degree-1] = nodes/2;
  
  int i = 0;
  while(1){
    int r = get_random((nodes-1)/2) + 1;
    if(used[r] != USED){
      adjacency[i] = r;
      adjacency[i+degree/2] = nodes - r;
      used[r] = USED;
      i++;
    }
    if(i == degree/2)
      break;
  }
}

static void init_used(const int nodes, bool used[nodes/2])
{
  for(int i=0;i<nodes/2;i++)
    used[i] = NOT_USED;

  used[0] = USED;
}

static void print_help(char *argv)
{
  ERROR("%s -N nodes -D degree [-o <output>] [-s <seed>] [-n <calcs>] [-w <max_temp>] [-c <min_temp>] [-A]\n", argv);
}

static void set_args(const int argc, char **argv, int *nodes, int *degree, char *outfname, bool *enable_output,
		     int *seed, long *ncalcs, double *max_temp, double *min_temp, bool *enable_ASPL_priority)
{
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
  int nodes = NOT_DEFINED, degree = NOT_DEFINED;
  int seed = 0, diameter, current_diameter, best_diameter, low_diameter;
  long sum, best_sum, ncalcs = 10000;
  double max_temp = 100, min_temp = 0.22, ASPL, current_ASPL, best_ASPL, low_ASPL;

  set_args(argc, argv, &nodes, &degree, outfname, &enable_output,
	   &seed, &ncalcs, &max_temp, &min_temp, &enable_ASPL_priority);
  if(nodes == NOT_DEFINED || degree == NOT_DEFINED)
    print_help(argv[0]);
  else if(nodes%2 == 1 && degree%2 == 1)
    ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);

  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  printf("Max, Min temperature = %f, %f\n", max_temp, min_temp);
  
  int lines           = (nodes * degree)/2;
  int (*edge)[2]      = malloc(sizeof(int)*lines*2);  // int edge[lines][2];
  int *adjacency      = malloc(sizeof(int)*degree);   // int adjacency[degree];
  int *best_adjacency = malloc(sizeof(int)*degree);   // int best_adjacency[degree];
  bool *used          = malloc(sizeof(bool)*nodes/2); // bool used[nodes/2];
  int symmetries      = nodes;

  ODP_Srand(seed);
  init_used(nodes, used);
  init_adjacency_general(nodes, degree, used, adjacency);

  ODP_Init_aspl_general_s(nodes, degree, NULL, symmetries); 
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);

  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL     = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int)*degree);

  ODP_Set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
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

      mutate_adjacency_general(nodes, degree, used, adjacency);
      ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
      
      if(diameter < best_diameter || (diameter == best_diameter && ASPL < best_ASPL)){
	best_diameter = diameter;
	best_sum      = sum;
	best_ASPL     = ASPL;
	memcpy(best_adjacency, adjacency, sizeof(int)*degree);
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
	restore_adjacency_general(nodes, degree, used, adjacency);
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
  printf("Time            = %f sec.\n", sa_time);
  printf("ASPL priority?  = %s\n", (enable_ASPL_priority)? "Yes" : "No");
  printf("Loop ?          = %s\n", (ODP_Check_loop(lines, edge))? "Yes" : "No");
  //  printf("Multiple Edges? = %s\n", (ODP_Check_multiple_edges(lines, edge))? "Yes" : "No");

  if(enable_output){
    ODP_Write_edge_general(lines, edge, outfname);
    printf("Generate ./%s\n", outfname);
  }
  
  free(edge);
  free(adjacency);
  free(best_adjacency);
  free(used);

  return 0;
}
