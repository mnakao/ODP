#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "apsp.h"

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  int diameter, low_diameter, width, height, rank;
  long sum;
  double ASPL, low_ASPL;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if(argc != 2){
    fprintf(stderr, "Please add an edge file. \"%s xxx.edges\"\n", argv[0]);
    exit(1);
  }
  
  int lines = apsp_get_lines(argv[1]);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  apsp_read_edge_grid(argv[1], &width, &height, edge);
  int nodes  = apsp_get_nodes(lines, edge);
  int degree = apsp_get_degree(nodes, lines, edge);

  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int adjacency[nodes][degree];
  apsp_set_adjacency(nodes, degree, lines, edge, adjacency);

  apsp_mpi_init(nodes, degree, NULL, MPI_COMM_WORLD);
  apsp_mpi_run(adjacency, &diameter, &sum, &ASPL);
  apsp_mpi_finalize();

  if(rank == 0){
    int length = apsp_get_length(lines, edge, height);
    apsp_set_lbounds_grid(width, height, degree, length, &low_diameter, &low_ASPL);
    printf("Width = %d, Height = %d, Length = %d, Degrees = %d\n", width, height, length, degree);
    printf("Diameter     = %d\n", diameter);
    printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
    printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
    printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);
  }
  
  free(edge);
  free(adjacency);
  MPI_Finalize();
  return 0;
}
