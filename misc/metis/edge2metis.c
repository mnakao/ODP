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
  int degree = ODP_Get_degree(nodes, lines, edge);
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree);
  ODP_Conv_edge2adjacency_general(nodes, lines, degree, edge, adjacency);
  
  printf("%d %d\n", nodes, (nodes*degree)/2);
  for(int i=0;i<nodes;i++){
    for(int j=0;j<degree;j++)
      printf("%d ", adjacency[i][j]+1);

    printf("\n");
  }
  
  return 0;
}
