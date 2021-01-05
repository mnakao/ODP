#include "common.h"

extern void apsp_cuda_init(const int nodes, const int degree, const int* num_degrees);
extern void apsp_cuda_run(const int* adjacency, int *diameter, long *sum, double *ASPL);
extern void apsp_cuda_finalize();

void apsp_all_cuda_run_general(const char *fname, int *nodes, int *degree,
			       int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
{
  int lines = apsp_get_lines(fname);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  apsp_read_edge_general(fname, edge);
  *nodes  = apsp_get_nodes(lines, edge);
  *degree = apsp_get_degree(*nodes, lines, edge);

  int (*adjacency)[*degree] = malloc(sizeof(int) * (*nodes) * (*degree)); // int adjacency[nodes][degree];
  apsp_conv_edge2adjacency(*nodes, lines, edge, (int *)adjacency);
  apsp_set_lbounds_general(*nodes, *degree, low_diameter, low_ASPL);

  int *num_degrees = malloc(sizeof(int) * (*nodes));
  apsp_set_degrees(*nodes, lines, edge, num_degrees);

  apsp_cuda_init(*nodes, *degree, num_degrees);
  apsp_cuda_run((int *)adjacency, diameter, sum, ASPL);
  apsp_cuda_finalize();

  free(num_degrees);
  free(adjacency);
  free(edge);
}

void apsp_all_cuda_run_grid(const char *fname, int *width, int *height, int *degree, int *length,
			    int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
{
  int lines = apsp_get_lines(fname);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  apsp_read_edge_grid(fname, width, height, edge);
  int nodes  = apsp_get_nodes(lines, edge);
  *degree = apsp_get_degree(nodes, lines, edge);

  int (*adjacency)[*degree] = malloc(sizeof(int) * nodes * (*degree)); // int adjacency[nodes][degree];
  apsp_conv_edge2adjacency(nodes, lines, edge, (int *)adjacency);
  *length = apsp_get_length(lines, edge, *height);
  apsp_set_lbounds_grid(*width, *height, *degree, *length, low_diameter, low_ASPL);

  int *num_degrees = malloc(sizeof(int) * nodes);
  apsp_set_degrees(nodes, lines, edge, num_degrees);
  
  apsp_cuda_init(nodes, *degree, num_degrees);
  apsp_cuda_run((int *)adjacency, diameter, sum, ASPL);
  apsp_cuda_finalize();

  free(num_degrees);
  free(adjacency);
  free(edge);
}
