#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "odp.h"

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
  
  int lines = ODP_Get_lines(argv[1]);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ODP_Read_edge_grid(argv[1], &width, &height, edge);
  int nodes  = ODP_Get_nodes(lines, edge);
  int degree = ODP_Get_degree(nodes, lines, edge);

  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int adjacency[nodes][degree];
  ODP_Conv_edge2adjacency_grid(width, height, lines, degree, edge, adjacency);

  ODP_Init_aspl_mpi_grid(width, height, degree, NULL, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // To measure time with ODP_PROFILE=1.
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  ODP_Finalize_aspl();

  if(rank == 0){
    int length = ODP_Get_length(lines, height, edge);
    ODP_Set_lbounds_grid(width, height, degree, length, &low_diameter, &low_ASPL);
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
