#include <stdio.h>
#include <stdlib.h>
#include "apsp.h"

int main(int argc, char *argv[])
{
  int width, height, degree, length, low_diameter, diameter;
  long sum;
  double low_ASPL, ASPL;

  if(argc != 2){
    fprintf(stderr, "Please add an edge file. \"%s xxx.edges\"\n", argv[0]);
    exit(1);
  }

  apsp_all_run_grid(argv[1], &width, &height, &degree, &length,
		    &low_diameter, &low_ASPL, &diameter, &sum, &ASPL);
    
  printf("Width = %d, Height = %d, Length = %d, Degrees = %d\n", width, height, length, degree);
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  int nodes = width * height;
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);

  return 0;
}
