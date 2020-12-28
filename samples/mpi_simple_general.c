#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "apsp.h"

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  int rank, nodes, degree, low_diameter, diameter;
  long sum;
  double ASPL, low_ASPL;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if(argc != 2){
    fprintf(stderr, "Please add an edge file. \"%s xxx.edges\"\n", argv[0]);
    exit(1);
  }
    
  apsp_all_mpi_run_general(argv[1], MPI_COMM_WORLD, &nodes, &degree,
			   &low_diameter, &low_ASPL, &diameter, &sum, &ASPL);

  if(rank == 0){
    apsp_set_lbounds_general(nodes, degree, &low_diameter, &low_ASPL);
    printf("Nodes = %d, Degrees = %d\n", nodes, degree);
    printf("Diameter     = %d\n", diameter);
    printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
    printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
    printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);
  }
  
  MPI_Finalize();
  return 0;
}
