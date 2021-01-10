#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "apsp.h"
#define MAX_FILENAME_LENGTH 256
#define ERROR(...) do{fprintf(stderr, __VA_ARGS__); exit(0);}while(0)

static void print_help(char *argv)
{
  ERROR("%s -N nodes -D degree [-o <output>] [-s <seed>] [-n <calcs>] [-w <max_temp>] [-c <min_temp>]\n", argv);
}

static void set_args(const int argc, char **argv, int *nodes, int *degree,
		     char *fname, int *seed, long *ncalcs, double *max_temp, double *min_temp)
{
  if(argc < 5)
    print_help(argv[0]);

  int result;
  while((result = getopt(argc,argv,"N:D:o:s:n:w:c:"))!=-1){
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
    default:
      print_help(argv[0]);
    }
  }
}

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

static bool accept(const int nodes, const int current_diameter, const int diameter,
		   const double current_ASPL, const double ASPL, const double temp)
{
  if(diameter < current_diameter){
    return true;
  }
  else if(diameter > current_diameter){
    return false;
  }
  else{ //  diameter == current_diameter
    if(ASPL <= current_ASPL){
      return true;
    }
    else{
      double diff = (current_ASPL-ASPL)*nodes*(nodes-1);
      if(exp(diff/temp) > uniform_rand()){
	return true;
      }
      else{
	return false;
      }
    }
  }
}

int main(int argc, char *argv[])
{
  int nodes, degree, seed = 0;
  double max_temp = 100, min_temp = 0.2;
  char *fname="general.edges";
  long sum, best_sum, ncalcs = 10000;
  int diameter, current_diameter, best_diameter, low_diameter;
  double ASPL, current_ASPL, best_ASPL, low_ASPL;

  set_args(argc, argv, &nodes, &degree, fname, &seed, &ncalcs, &max_temp, &min_temp);
  if(nodes%2 == 1 && degree%2 == 1)
    ERROR("Invalid nodes(%d) or degree(%d)\n", nodes, degree);
  
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  printf("Max, Min temperature = %f, %f\n", max_temp, min_temp);
  
  int lines = (nodes * degree)/2;
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int adjacency[nodes][degree];
  int (*best_adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int best_adjacency[nodes][degree];
  
  apsp_srand(seed);
  apsp_generate_random_general(nodes, degree, edge);
  apsp_conv_edge2adjacency(nodes, lines, edge, adjacency);
  
  apsp_run_init(nodes, degree, NULL);
  apsp_run(adjacency, &diameter, &sum, &ASPL);
  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL     = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int)*nodes*degree);

  apsp_set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
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
      
      apsp_mutate_adjacency_general(nodes, degree, NULL, adjacency);
      apsp_run(adjacency, &diameter, &sum, &ASPL);
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
      
      if(accept(nodes, current_diameter, diameter, current_ASPL, ASPL, temp)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
      }
      else{
	apsp_restore_adjacency(adjacency);
      }
      temp *= cooling_rate;
    }
  }
  
  apsp_run_finalize();
  apsp_conv_adjacency2edge(nodes, degree, NULL, best_adjacency, edge);
  printf("---\n");
  printf("Diameter       = %d\n", best_diameter);
  printf("Diameter Gap   = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL           = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap       = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("Multiple edges = %s\n", (apsp_check_multiple_edges(lines, edge))? "Exist" : "None");
  printf("Loop           = %s\n", (apsp_check_loop(lines, edge))? "Exist" : "None");

  apsp_write_edge_general(lines, edge, fname);
  printf("Generate ./%s\n", fname);
  
  free(edge);
  free(adjacency);
  free(best_adjacency);

  return 0;
}
