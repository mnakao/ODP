#include "common.h"

extern void ODP_Init_aspl_general(const int nodes, const int degree, const int* num_degrees);
extern void ODP_Init_aspl_grid(const int width, const int height, const int degree, const int* num_degrees);
extern void ODP_Set_aspl(const int* restrict adjacency, int *diameter, long *sum, double *ASPL);
extern void ODP_Finalize_aspl();

void ODP_Set_aspl_general(const char *fname, int *nodes, int *degree,
                          int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
{
  int lines = ODP_Get_lines(fname);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ODP_Read_edge_general(fname, edge);
  *nodes  = ODP_Get_nodes(lines, edge);
  *degree = ODP_Get_degree(*nodes, lines, edge);

  int (*adjacency)[*degree] = malloc(sizeof(int) * (*nodes) * (*degree)); // int adjacency[nodes][degree];
  ODP_Conv_edge2adjacency_general(*nodes, lines, *degree, edge, (int *)adjacency);
  ODP_Set_lbounds_general(*nodes, *degree, low_diameter, low_ASPL);

  int *num_degrees = malloc(sizeof(int) * (*nodes));
  ODP_Set_degrees(*nodes, lines, edge, num_degrees);
  
  ODP_Init_aspl_general(*nodes, *degree, num_degrees);
  ODP_Set_aspl((int*)adjacency, diameter, sum, ASPL);
  ODP_Finalize_aspl();

  free(num_degrees);
  free(adjacency);
  free(edge);
}

void ODP_Set_aspl_grid(const char *fname, int *width, int *height, int *degree, int *length,
                       int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
{
  int lines = ODP_Get_lines(fname);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ODP_Read_edge_grid(fname, width, height, edge);
  int nodes = (*width) * (*height);
  *degree = ODP_Get_degree(nodes, lines, edge);

  int (*adjacency)[*degree] = malloc(sizeof(int) * nodes * (*degree)); // int adjacency[nodes][degree];
  ODP_Conv_edge2adjacency_grid(*width, *height, lines, *degree, edge, (int *)adjacency);
  *length = ODP_Get_length(lines, edge, *height);
  ODP_Set_lbounds_grid(*width, *height, *degree, *length, low_diameter, low_ASPL);

  int *num_degrees = malloc(sizeof(int) * nodes);
  ODP_Set_degrees(nodes, lines, edge, num_degrees);

  ODP_Init_aspl_grid(*width, *height, *degree, num_degrees);
  ODP_Set_aspl((int*)adjacency, diameter, sum, ASPL);
  ODP_Finalize_aspl();

  free(num_degrees);
  free(adjacency);
  free(edge);
}
