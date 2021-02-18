#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "odp.h"
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)

static void set_args(const int argc, char **argv, int *nodes, int *degree, int *groups, int *seed)
{
  int result;
  while((result = getopt(argc,argv,"n:d:g:s:"))!=-1){
    switch(result){
    case 'n':
      *nodes = atoi(optarg);
      if(*nodes <= 0)
        ERROR("-n value > 0\n");
      break;
    case 'd':
      *degree = atoi(optarg);
      if(*degree <= 0)
        ERROR("-d value > 0\n");
      break;
    case 'g':
      *groups = atoi(optarg);
      if(*groups <= 0)
        ERROR("-g value > 0\n");
      break;
    case 's':
      *seed = atoi(optarg);
      if(*seed < 0)
        ERROR("-s value >= 0\n");
      break;
    }
  }
}

int main(int argc, char *argv[])
{
  int based_nodes = 10, degree = 3, groups = 2, seed = 1;
  set_args(argc, argv, &based_nodes, &degree, &groups, &seed);
  int nodes = based_nodes * groups;
  int lines = (nodes * degree) / 2;
  int edge[lines][2];

  //  if(nodes%2 != 0 && degree%2 != 0)
  //    ERROR("Error! 123\n");

  ODP_Srand(seed);
  ODP_Generate_random_general_s(nodes, degree, groups, edge);
  ODP_Print_edge_general(lines, edge);
  
  return 0;
}
