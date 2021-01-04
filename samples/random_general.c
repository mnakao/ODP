#include <stdio.h>
#include <stdlib.h>
#include "apsp.h"

int main()
{
  int nodes = 100, degree = 5, seed = 1, diameter, low_diameter;
  int lines = (nodes * degree)/2;
  long sum;
  double ASPL, low_ASPL;
  char *fname="general.edges";

  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  apsp_random_general(nodes, degree, seed, edge);
  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int adjacency[nodes][degree];
  apsp_set_adjacency(nodes, degree, lines, edge, adjacency);

  apsp_init(nodes, degree, NULL);
  apsp_run(adjacency, &diameter, &sum, &ASPL);
  apsp_finalize();

  apsp_set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);

  apsp_write_edge_general(lines, edge, fname);
  printf("Generate ./%s\n", fname);
  
  free(edge);
  free(adjacency);

  return 0;
}
