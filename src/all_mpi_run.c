#include "common.h"
#include <mpi.h>
extern void apsp_mpi_init(const int nodes, const int degree, const int* restrict num_degrees, MPI_Comm comm);
extern void apsp_mpi_run(const int* restrict adjacency, int *diameter, long *sum, double *ASPL);
extern void apsp_mpi_finalize();

void apsp_all_mpi_run_general(const char *fname, const MPI_Comm comm, int *nodes, int *degree,
			      int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
{
  int lines = apsp_get_lines(fname);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  apsp_read_edge_general(fname, edge);
  *nodes  = apsp_get_nodes(lines, edge);
  *degree = apsp_get_degree(*nodes, lines, edge);

  int (*adjacency)[*degree] = malloc(sizeof(int) * (*nodes) * (*degree)); // int adjacency[nodes][degree];
  apsp_conv_edge2adjacency(*nodes, *degree, lines, edge, (int *)adjacency);
  apsp_set_lbounds_general(*nodes, *degree, low_diameter, low_ASPL);

  int *num_degrees = malloc(sizeof(int) * (*nodes));
  apsp_set_degrees(*nodes, lines, edge, num_degrees);

  apsp_mpi_init(*nodes, *degree, num_degrees, comm);
  apsp_mpi_run((int*)adjacency, diameter, sum, ASPL);
  apsp_mpi_finalize();

  free(num_degrees);
  free(adjacency);
  free(edge);
}

void apsp_all_mpi_run_grid(const char *fname, const MPI_Comm comm, int *width, int *height, int *degree, int *length,
			   int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
{
  int lines = apsp_get_lines(fname);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  apsp_read_edge_grid(fname, width, height, edge);
  int nodes  = apsp_get_nodes(lines, edge);
  *degree = apsp_get_degree(nodes, lines, edge);

  int (*adjacency)[*degree] = malloc(sizeof(int) * nodes * (*degree)); // int adjacency[nodes][degree];
  apsp_conv_edge2adjacency(nodes, *degree, lines, edge, (int *)adjacency);
  *length = apsp_get_length(lines, edge, *height);
  apsp_set_lbounds_grid(*width, *height, *degree, *length, low_diameter, low_ASPL);

  int *num_degrees = malloc(sizeof(int) * nodes);
  apsp_set_degrees(nodes, lines, edge, num_degrees);

  apsp_mpi_init(nodes, *degree, num_degrees, comm);
  apsp_mpi_run((int*)adjacency, diameter, sum, ASPL);
  apsp_mpi_finalize();
  
  free(num_degrees);
  free(adjacency);
  free(edge);
}
