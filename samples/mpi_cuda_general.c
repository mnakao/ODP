#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "odp.h"

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  int diameter, low_diameter, rank;
  long sum;
  double ASPL, low_ASPL;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if(argc != 2){
    fprintf(stderr, "Please add an edge file. \"%s xxx.edges\"\n", argv[0]);
    exit(1);
  }
    
  int lines = ODP_Get_lines(argv[1]);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ODP_Read_edge_general(argv[1], edge);
  int nodes  = ODP_Get_nodes(lines, edge);
  int degree = ODP_Get_degree(nodes, lines, edge);

  int (*adjacency)[degree] = malloc(sizeof(int) * nodes * degree); // int adjacency[nodes][degree];
  ODP_Conv_edge2adjacency(nodes, lines, degree, edge, adjacency);
  
  ODP_Init_aspl_mpi_cuda_general(nodes, degree, NULL, MPI_COMM_WORLD);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  ODP_Finalize_aspl();

  if(rank == 0){
    ODP_Set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
    printf("Nodes = %d, Degrees = %d\n", nodes, degree);
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
