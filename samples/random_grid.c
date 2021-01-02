#include <stdio.h>
#include <stdlib.h>
#include "apsp.h"

int main()
{
  int width = 12, height = 8, degree = 3, length = 3, seed = 1;
  int nodes = width * height;
  int lines = (nodes*degree)/2, diameter, low_diameter;
  long sum;
  double ASPL, low_ASPL;
  char *fname="grid.edges";

  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  apsp_random_grid(width, height, degree, length, seed, edge);
  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int adjacency[nodes][degree];
  apsp_set_adjacency(nodes, degree, lines, edge, adjacency);

  apsp_init(nodes, degree, NULL);
  apsp_run(adjacency, &diameter, &sum, &ASPL);
  apsp_finalize();

  apsp_set_lbounds_grid(width, height, degree, length, &low_diameter, &low_ASPL);
  printf("Width = %d, Height = %d, Length = %d, Degrees = %d\n", width, height, length, degree);
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);

  apsp_output_edge_grid(fname, lines, height, edge);
  printf("Generate ./%s\n", fname);
  
  free(edge);
  free(adjacency);

  return 0;
}
