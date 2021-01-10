#ifndef APSP_H_INCLUDED
#define APSP_H_INCLUDED
#include <stdbool.h>

extern void apsp_run_init(int nodes, int degree, int* num_degrees);
extern void apsp_run_init_s(int nodes, int degree, int* num_degrees, int symmetries);
extern void apsp_run_finalize();
extern void apsp_run(void *adjacency, int *diameter, long *sum, double *ASPL);
extern void apsp_all_run_general(char *fname, int *nodes, int *degree,
				 int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);
extern void apsp_all_run_grid(char *fname, int *width, int *height, int *degree, int *length,
			      int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);

extern void apsp_cuda_run_init(int nodes, int degree, int* num_degrees);
extern void apsp_cuda_run_init_s(int nodes, int degree, int* num_degrees, int symmetries);
extern void apsp_cuda_run_finalize();
extern void apsp_cuda_run(void *adjacency, int *diameter, long *sum, double *ASPL);
extern void apsp_all_cuda_run_general(char *fname, int *nodes, int *degree, int *low_diameter,
				      double *low_ASPL, int *diameter, long *sum, double *ASPL);
extern void apsp_all_cuda_run_grid(char *fname, int *width, int *height, int *degree, int *length,
				   int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);

#ifdef MPI_VERSION
#include <mpi.h>
extern void apsp_mpi_run_init(int nodes, int degree, int* num_degrees, MPI_Comm comm);
extern void apsp_mpi_run_init_s(int nodes, int degree, int* num_degrees, MPI_Comm comm, int symmetries);
extern void apsp_mpi_run_finalize();
extern void apsp_mpi_run(void *adjacency, int *diameter, long *sum, double *ASPL);
extern void apsp_all_mpi_run_general(char *fname, MPI_Comm comm, int *nodes, int *degree,
				     int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);
extern void apsp_all_mpi_run_grid(char *fname, MPI_Comm comm, int *width, int *height, int *degree, int *length,
				  int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);

extern void apsp_mpi_cuda_run_init(int nodes, int degree, int* num_degrees, MPI_Comm comm);
extern void apsp_mpi_cuda_run_init_s(int nodes, int degree, int* num_degrees, MPI_Comm comm, int symmetries);
extern void apsp_mpi_cuda_run_finalize();
extern void apsp_mpi_cuda_run(void *adjacency, int *diameter, long *sum, double *ASPL);
extern void apsp_all_mpi_cuda_run_general(char *fname, MPI_Comm comm, int *nodes, int *degree,
					  int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);
extern void apsp_all_mpi_cuda_run_grid(char *fname, MPI_Comm comm, int *width, int *height, int *degree, int *length,
				       int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);
#endif

// APSP utilities
extern void apsp_set_lbounds_general(int nodes, int degree, int *low_diameter, double *low_ASPL);
extern void apsp_set_lbounds_grid(int m, int n, int degree, int length, int *low_diameter, double *low_ASPL);
extern void apsp_set_degrees(int nodes, int lines, int (*edge)[2], int* num_degrees);
extern int  apsp_get_lines(char *fname);
extern int  apsp_get_nodes(int lines, int (*edge)[2]);
extern int  apsp_get_degree(int nodes, int lines, int (*edge)[2]);
extern int  apsp_get_length(int lines, int (*edge)[2], int height);
extern bool apsp_check_general(char *fname);
extern bool apsp_check_multiple_edges(int lines, int (*edge)[2]);
extern bool apsp_check_loop(int lines, int (*edge)[2]);
extern void apsp_read_edge_general(char* fname, int (*edge)[2]);
extern void apsp_read_edge_grid(char *fname, int *w, int *h, int (*edge)[2]);
extern void apsp_write_edge_general(int lines, int (*edge)[2], char *fname);
extern void apsp_write_edge_grid(int lines, int height, int (*edge)[2], char *fname);
extern void apsp_conv_edge2adjacency(int nodes, int lines, int (*edge)[2], void *adjacency);
extern void apsp_conv_adjacency2edge(int nodes, int degree, int *num_degrees, void *adjacency, int (*edge)[2]);
extern void apsp_conv_edge2adjacency_s(int nodes, int lines, int (*edge)[2], int symmetries, void *adjacency);
extern void apsp_generate_random_general(int nodes, int degree, int (*edge)[2]);
extern void apsp_generate_random_grid(int width, int height, int degree, int length, int (*edge)[2]);
extern void apsp_srand(unsigned int seed);
extern void apsp_mutate_adjacency_general(int nodes, int degree, int *num_degrees, void *adjacency);
extern void apsp_restore_adjacency(void *adjacency);
extern void apsp_print_adjacency(int nodes, int degree, int *num_degrees, void *adjacency);
extern void apsp_print_edge(int lines, int (*edge)[2]);
#endif
