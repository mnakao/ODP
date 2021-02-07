#include <stdio.h>
#include <stdlib.h>
#include "odp.h"

int main()
{
  int width = 12, height = 8, degree = 3, length = 3;
  int nodes = width * height;
  int lines = (nodes*degree)/2, diameter, low_diameter;
  unsigned int seed = 1;
  long sum;
  double ASPL, low_ASPL;
  char *fname="grid.edges";

  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ODP_Srand(seed);
  ODP_Generate_random_grid(width, height, degree, length, edge);
  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int adjacency[nodes][degree];
  ODP_Conv_edge2adjacency_grid(width, height, lines, degree, edge, adjacency);

  ODP_Init_aspl_grid(width, height, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  ODP_Finalize_aspl();

  ODP_Set_lbounds_grid(width, height, degree, length, &low_diameter, &low_ASPL);
  printf("Width = %d, Height = %d, Length = %d, Degrees = %d\n", width, height, length, degree);
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);

  ODP_Write_edge_grid(lines, height, edge, fname);
  printf("Generate ./%s\n", fname);
  
  free(edge);
  free(adjacency);

  return 0;
}
