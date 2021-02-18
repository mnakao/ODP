#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "odp.h"
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)

static void set_args(const int argc, char **argv, int *width, int *height, int *length, int *degree, int *groups, int *seed)
{
  int result;
  while((result = getopt(argc,argv,"w:h:l:d:g:s:"))!=-1){
    switch(result){
    case 'w':
      *width = atoi(optarg);
      if(*width <= 0)
        ERROR("-w value > 0\n");
      break;
    case 'h':
      *height = atoi(optarg);
      if(*height <= 0)
        ERROR("-h value > 0\n");
      break;
    case 'l':
      *length = atoi(optarg);
      if(*length <= 0)
        ERROR("-l value > 0\n");
      break;
    case 'd':
      *degree = atoi(optarg);
      if(*degree <= 0)
        ERROR("-d value > 0\n");
      break;
    case 'g':
      *groups = atoi(optarg);
      if(*groups <= 0)
        ERROR("-g value > 0\n");
      break;
    case 's':
      *seed = atoi(optarg);
      if(*seed < 0)
        ERROR("-s value >= 0\n");
      break;
    }
  }
}

int main(int argc, char *argv[])
{
  int width = 6, height = 6, degree = 3, length = 3, seed = 1, symmetries = 2;
  set_args(argc, argv, &width, &height, &length, &degree, &symmetries, &seed);
  int nodes = width * height;
  int lines = nodes*degree/2;
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];

  ODP_Srand(seed);
  ODP_Generate_random_grid_s(width, height, degree, length, symmetries, edge);
  ODP_Print_edge_grid(lines, height, edge);
  //  printf("%d\n", ODP_Get_length(lines, height, edge));
  //  ODP_Generate_random_general_s(nodes, degree, seed, symmetries, edge);
  //  ODP_Print_edge_general(lines, edge);

  //  int based_nodes = nodes/symmetries;
  //  int (*adjacency)[degree] = malloc(sizeof(int) * based_nodes * degree); // int adjacency[based_nodes][degree];
  //  int edge2[50][2] = {{0, 2},{0, 18},{1, 3},{1, 4},{1, 8},{2, 5},{2, 8},{3, 9},{3, 10},{4, 10},{5, 17},{6, 24},{6, 7},{7, 14},{7, 12},{8, 26},{9, 27},{9, 20},{10, 15},{11, 23},{11, 29},{11, 16},{12, 13},{12, 24},{13, 16},{13, 21},{14, 21},{14, 22},{15, 20},{15, 26},{16, 17},
		      //		      {17, 35},{18, 19},{18, 30},{19, 24},{19, 22},{20, 25},{21, 28},{22, 23},{23, 28},{25, 31},{25, 32},{26, 32},{27, 33},{27, 34},{28, 29},{30, 33},{31, 34},{32, 34},{33 ,35}};
  //  ODP_Conv_edge2adjacency_grid_s(width, height, lines-4, degree, edge2, symmetries, adjacency);
  
  //  ODP_Conv_adjacency2edge_grid_s(width, height, degree, NULL, adjacency, symmetries, edge);
  //  ODP_Print_edge(lines, edge);
  //  ODP_Print_adjacency(based_nodes, degree, NULL, adjacency);
  /*
  ODP_Init_aspl(nodes, degree, NULL);
  ODP_Set_aspl(adjacency, &diameter, &sum, &ASPL);
  ODP_Finalize_aspl();

  int length = ODP_Get_length(lines, height, edge);
  ODP_Set_lbounds_grid(width, height, degree, length, &low_diameter, &low_ASPL);
  printf("Width = %d, Height = %d, Length = %d, Degrees = %d\n", width, height, length, degree);
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);
  */
  return 0;
}
