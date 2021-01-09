#include <stdio.h>
#include "apsp.h"

int main()
{
  int nodes = 12, degree = 3, lines = 15, diameter, low_diameter;
  long sum;
  double ASPL, low_ASPL;
  int edge[][2] = {{0,10},{0,3},{0,4},{1,8},{1,3},{1,7},{2,8},{2,9},{2,6},{3,5},{4,9},{4,10},{5,11},{5,6},{6,7}};
  
  int adjacency[nodes][degree], num_degrees[nodes];
  apsp_conv_edge2adjacency(nodes, lines, edge, adjacency);
  // adjacency[][] = {{10,3,4}, {8,3,7}, {8,9,6}, {0,1,5}, {0,9,10}, {3,11,6}, {2,5,7}, {1,6}, {1,2}, {2,4}, {0,4}, {5}};
  apsp_set_degrees(nodes, lines, edge, num_degrees); // If it is a non-regular graph, count the number of edges that each vertex has.
  // num_degrees[] = {3,3,3, 3,3,3, 3,2,2, 2,2,1};

  // In a regular graph, the third argument is NULL.
  // In a non-regular graph, the third argument is the number of edges that each vertex has.
  apsp_run_init(nodes, degree, num_degrees);
  apsp_run(adjacency, &diameter, &sum, &ASPL);
  apsp_run_finalize();
  
  apsp_set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);

  return 0;
}
