#include <stdio.h>
#include <stdlib.h>
#include "odp.h"

int main()
{
  int nodes = 100, degree = 5, diameter, low_diameter;
  int lines = (nodes * degree)/2;
  unsigned int seed = 1;
  long sum;
  double ASPL, low_ASPL;
  char *fname="general.edges";

  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ODP_Srand(seed);
  ODP_Generate_random_general(nodes, degree, edge);
  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int adjacency[nodes][degree];
  ODP_Conv_edge2adjacency(nodes, lines, degree, edge, adjacency);

  ODP_Init_aspl(nodes, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  ODP_Finalize_aspl();

  ODP_Set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);

  ODP_Write_edge_general(lines, edge, fname);
  printf("Generate ./%s\n", fname);
  
  free(edge);
  free(adjacency);

  return 0;
}
