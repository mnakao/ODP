#include <stdio.h>
#include <stdlib.h>
#include "odp.h"

int main(int argc, char *argv[])
{
  if(argc != 2){
    fprintf(stderr, "%s edge_file\n", argv[0]);
    exit(1);
  }
  char *infname  = argv[1];
  int lines      = ODP_Get_lines(infname);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ODP_Read_edge_general(infname, edge);
  int nodes  = ODP_Get_nodes(lines, edge);
  
  printf("%d %d\n", lines, nodes);
  for(int i=0;i<lines;i++)
    printf("%d %d\n", edge[i][0]+1, edge[i][1]+1);

  return 0;
}
