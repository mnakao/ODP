#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include "odp.h"
#define NOT_DEFINED -1
#define MAX_FILENAME_LENGTH 256
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)

static void print_help(char *argv)
{
  ERROR("%s -H hosts -R radix -f input\n", argv);
}

static void set_args(const int argc, char **argv, int *hosts, int *radix, char **infname)
{
  int result;
  while((result = getopt(argc,argv,"H:R:f:"))!=-1){
    switch(result){
    case 'H':
      *hosts = atoi(optarg);
      if(*hosts <= 0)
        ERROR("-H value > 0\n");
      break;
    case 'R':
      *radix = atoi(optarg);
      if(*radix <= 0)
        ERROR("-R value > 0\n");
      break;
    case 'f':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Input filename is long (%s).\n", optarg);
      *infname = malloc(MAX_FILENAME_LENGTH);
      strcpy(*infname, optarg);
      break;
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  int hosts = NOT_DEFINED, radix = NOT_DEFINED;
  char *infname = NULL;
  
  set_args(argc, argv, &hosts, &radix, &infname);
  int lines = ODP_Get_lines(infname);
  int (*odp_edge)[2] = malloc(sizeof(int)*lines*2);
  ODP_Read_edge_general(infname, odp_edge);
  int nodes  = ODP_Get_nodes(lines, odp_edge);
  int degree = ODP_Get_degree(nodes, lines, odp_edge);
  int switches = nodes;
  if(hosts == NOT_DEFINED || radix == NOT_DEFINED) print_help(argv[0]);
  else if(hosts%switches != 0) ERROR("hosts%switches != 0\n");
  else if(degree >= radix) ERROR("degree >= radix\n");
  
  int (*orp_edge)[2] = malloc(sizeof(int)*(hosts*2));
  int k = 0;
  for(int i=0;i<switches;i++){
    for(int j=0;j<hosts/switches;j++){
      orp_edge[k][0] = k;
      orp_edge[k][1] = hosts + i;
      k++;
    }
  }
  
  printf("%d %d %d\n", hosts, switches, radix);
  for(int i=0;i<hosts;i++)
    printf("%d %d\n", orp_edge[i][0], orp_edge[i][1]);
  for(int i=0;i<lines;i++)
    printf("%d %d\n", odp_edge[i][0]+hosts, odp_edge[i][1]+hosts);
  
  return 0;
}
