#include <stdio.h>
#include "odp.h"

int main()
{
  int nodes = 12, degree = 3, lines = 18, diameter, low_diameter;
  long sum;
  double ASPL, low_ASPL;
  int edge[][2] = {{0,10},{0,3},{0,4},{1,8},{1,3},{1,7},{2,8},{2,9},{2,6},{3,5},{4,9},{4,10},{5,11},{5,6},{6,7},{7,11},{8,11},{9,10}};
  // ../data/general/n12d3.random.edges
  
  int adjacency[nodes][degree];
  ODP_Conv_edge2adjacency(nodes, lines, edge, adjacency);
  // adjacency[][] = {{10,3,4}, {8,3,7}, {8,9,6}, {0,1,5}, {0,9,10}, {3,11,6}, {2,5,7}, {1,6,11}, {1,2,11}, {2,4,10}, {0,4,9}, {5,7,8}};

  ODP_Init_aspl(nodes, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  ODP_Finalize_aspl();

  ODP_Set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);

  return 0;
}
