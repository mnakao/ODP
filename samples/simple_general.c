#include <stdio.h>
#include <stdlib.h>
#include "odp.h"

int main(int argc, char *argv[])
{
  int nodes, degree, low_diameterm, diameter, low_diameter;
  long sum;
  double low_ASPL, ASPL;

  if(argc != 2){
    fprintf(stderr, "Please add an edge file. \"%s xxx.edges\"\n", argv[0]);
    exit(1);
  }
  
  ODP_Set_aspl_general(argv[1], &nodes, &degree, &low_diameter, &low_ASPL, &diameter, &sum, &ASPL);
   
  printf("Nodes = %d, Degrees = %d\n", nodes, degree);
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);
  
  return 0;
}
