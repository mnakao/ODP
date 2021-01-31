#ifndef APSP_H_INCLUDED
#define APSP_H_INCLUDED
#include <stdbool.h>

extern void ODP_Init_aspl(int nodes, int degree, int* num_degrees);
extern void ODP_Init_aspl_s(int nodes, int degree, int* num_degrees, int symmetries);
extern void ODP_Finalize_aspl();
extern void ODP_Set_aspl(void *adjacency, int *diameter, long *sum, double *ASPL);
extern void ODP_Set_aspl_general(char *fname, int *nodes, int *degree,
				 int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);
extern void ODP_Set_aspl_grid(char *fname, int *width, int *height, int *degree, int *length,
			      int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);

extern void ODP_Init_aspl_cuda(int nodes, int degree, int* num_degrees);
extern void ODP_Init_aspl_cuda_s(int nodes, int degree, int* num_degrees, int symmetries);
extern void ODP_Finalize_aspl_cuda();
extern void ODP_Set_aspl_cuda(void *adjacency, int *diameter, long *sum, double *ASPL);
extern void ODP_Set_aspl_cuda_general(char *fname, int *nodes, int *degree, int *low_diameter,
				      double *low_ASPL, int *diameter, long *sum, double *ASPL);
extern void ODP_Set_aspl_cuda_grid(char *fname, int *width, int *height, int *degree, int *length,
				   int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);

#ifdef MPI_VERSION
#include <mpi.h>
extern void ODP_Init_aspl_mpi(int nodes, int degree, int* num_degrees, MPI_Comm comm);
extern void ODP_Init_aspl_mpi_s(int nodes, int degree, int* num_degrees, MPI_Comm comm, int symmetries);
extern void ODP_Finalize_aspl_mpi();
extern void ODP_Set_aspl_mpi(void *adjacency, int *diameter, long *sum, double *ASPL);
extern void ODP_Set_aspl_mpi_general(char *fname, MPI_Comm comm, int *nodes, int *degree,
				     int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);
extern void ODP_Set_aspl_mpi_grid(char *fname, MPI_Comm comm, int *width, int *height, int *degree, int *length,
				  int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);

extern void ODP_Init_aspl_mpi_cuda(int nodes, int degree, int* num_degrees, MPI_Comm comm);
extern void ODP_Init_aspl_mpi_cuda_s(int nodes, int degree, int* num_degrees, MPI_Comm comm, int symmetries);
extern void ODP_Finalize_aspl_mpi_cuda();
extern void ODP_Set_aspl_mpi_cuda(void *adjacency, int *diameter, long *sum, double *ASPL);
extern void ODP_Set_aspl_mpi_cuda_general(char *fname, MPI_Comm comm, int *nodes, int *degree,
					  int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);
extern void ODP_Set_aspl_mpi_cuda_grid(char *fname, MPI_Comm comm, int *width, int *height, int *degree, int *length,
				       int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL);
#endif

// APSP utilities
extern void ODP_Set_lbounds_general(int nodes, int degree, int *low_diameter, double *low_ASPL);
extern void ODP_Set_lbounds_grid(int m, int n, int degree, int length, int *low_diameter, double *low_ASPL);
extern void ODP_Set_degrees(int nodes, int lines, int (*edge)[2], int* num_degrees);
extern int  ODP_Get_lines(char *fname);
extern int  ODP_Get_nodes(int lines, int (*edge)[2]);
extern int  ODP_Get_degree(int nodes, int lines, int (*edge)[2]);
extern int  ODP_Get_length(int lines, int (*edge)[2], int height);
extern bool ODP_Check_general(char *fname);
extern bool ODP_Check_multiple_edges(int lines, int (*edge)[2]);
extern bool ODP_Check_loop(int lines, int (*edge)[2]);
extern void ODP_Read_edge_general(char* fname, int (*edge)[2]);
extern void ODP_Read_edge_grid(char *fname, int *w, int *h, int (*edge)[2]);
extern void ODP_Write_edge_general(int lines, int (*edge)[2], char *fname);
extern void ODP_Write_edge_grid(int lines, int height, int (*edge)[2], char *fname);
extern void ODP_Conv_edge2adjacency(int nodes, int lines, int degree, int (*edge)[2], void *adjacency);
extern void ODP_Conv_edge2adjacency_general_s(int nodes, int lines, int degree, int (*edge)[2], int symmetries, void *adjacency);
extern void ODP_Conv_edge2adjacency_grid_s(int width, int height, int lines, int degree, int (*edge)[2], int symmetries, int (*adjacency)[degree]);
extern void ODP_Conv_adjacency2edge(int nodes, int degree, int *num_degrees, void *adjacency, int (*edge)[2]);
extern void ODP_Conv_adjacency2edge_general_s(int nodes, int degree,int *num_degrees, void *adjacency, int symmetries, int (*edge)[2]);
extern void ODP_Conv_adjacency2edge_grid_s(int width, int height, int degree,int *num_degrees, void *adjacency, int symmetries, int (*edge)[2]);
extern void ODP_Generate_random_general(int nodes, int degree, unsigned int seed, int (*edge)[2]);
extern void ODP_Generate_random_general_s(int nodes, int degree, unsigned int seed, int symmetries, int (*edge)[2]);
extern void ODP_Generate_random_grid(int width, int height, int degree, int length, unsigned int seed, int (*edge)[2]);
extern void ODP_Generate_random_grid_s(int width, int height, int degree, int length, unsigned int seed, int symmetries, int (*edge)[2]);
extern void ODP_Print_adjacency(int nodes, int degree, int *num_degrees, void *adjacency);
extern void ODP_Print_edge_general(int lines, int (*edge)[2]);
extern void ODP_Print_edge_grid(int lines, int height, int (*edge)[2]);
#endif
