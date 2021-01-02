#ifndef APSP_H_INCLUDED
#define APSP_H_INCLUDED
#include <stdbool.h>

extern void apsp_init(int nodes, int degree, int* num_degrees);
extern void apsp_init_s(int nodes, int degree, int* num_degrees, int groups);
extern void apsp_finalize();
extern void apsp_run(void* adjacency, int *diameter, long *sum, double *ASPL);
extern void apsp_all_run_general(char *fname, int *nodes, int *degree,
				 int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);
extern void apsp_all_run_grid(char *fname, int *width, int *height, int *degree, int *length,
			      int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);

extern void apsp_cuda_init(int nodes, int degree, int* num_degrees);
extern void apsp_cuda_init_s(int nodes, int degree, int* num_degrees, int groups);
extern void apsp_cuda_finalize();
extern void apsp_cuda_run(void* adjacency, int *diameter, long *sum, double *ASPL);
extern void apsp_all_cuda_run_general(char *fname, int *nodes, int *degree, int *low_diameter,
				      double *low_ASPL, int *diameter, long *sum, double *ASPL);
extern void apsp_all_cuda_run_grid(char *fname, int *width, int *height, int *degree, int *length,
				   int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);

#ifdef MPI_VERSION
#include <mpi.h>
extern void apsp_mpi_init(int nodes, int degree, int* num_degrees, MPI_Comm comm);
extern void apsp_mpi_init_s(int nodes, int degree, int* num_degrees, MPI_Comm comm, int groups);
extern void apsp_mpi_finalize();
extern void apsp_mpi_run(void* adjacency, int *diameter, long *sum, double *ASPL);
extern void apsp_all_mpi_run_general(char *fname, MPI_Comm comm, int *nodes, int *degree,
				     int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);
extern void apsp_all_mpi_run_grid(char *fname, MPI_Comm comm, int *width, int *height, int *degree, int *length,
				  int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);

extern void apsp_mpi_cuda_init(int nodes, int degree, int* num_degrees, MPI_Comm comm);
extern void apsp_mpi_cuda_init_s(int nodes, int degree, int* num_degrees, MPI_Comm comm, int groups);
extern void apsp_mpi_cuda_finalize();
extern void apsp_mpi_cuda_run(void* adjacency, int *diameter, long *sum, double *ASPL);
extern void apsp_all_mpi_cuda_run_general(char *fname, MPI_Comm comm, int *nodes, int *degree,
					  int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);
extern void apsp_all_mpi_cuda_run_grid(char *fname, MPI_Comm comm, int *width, int *height, int *degree, int *length,
				       int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);
#endif

// APSP utilities
extern void apsp_set_lbounds_general(int nodes, int degree, int *low_diameter, double *low_ASPL);
extern void apsp_set_lbounds_grid(int m, int n, int degree, int length, int *low_diameter, double *low_ASPL);
extern void apsp_set_edge_general(char* fname, int (*edge)[2]);
extern void apsp_set_edge_grid(char *fname, int *w, int *h, int (*edge)[2]);
extern void apsp_set_adjacency(int nodes, int degree, int lines, int edge[lines][2], int adjacency[nodes][degree]);
extern void apsp_set_degrees(int nodes, int lines, int edge[lines][2], int* num_degrees);
extern int  apsp_get_lines(char* fname);
extern int  apsp_get_nodes(int lines, int edge[lines][2]);
extern int  apsp_get_degree(int nodes, int lines, int edge[lines][2]);
extern int  apsp_get_length(int lines, int edge[lines][2], int height);
extern bool apsp_check_general(char *fname);
extern bool apsp_check_duplicated_edge(int lines, int edge[lines][2]);
extern bool apsp_check_loop(int lines, int edge[lines][2]);
extern void apsp_random_general(int nodes, int degree, unsigned int seed, void *edge);
extern void apsp_output_edge_general(char *fname, int lines, int (*edge)[2]);
#endif
