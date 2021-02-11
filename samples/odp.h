#ifndef APSP_H_INCLUDED
#define APSP_H_INCLUDED
#include <stdbool.h>

extern void ODP_Init_aspl_general(int nodes, int degree, int* num_degrees);
extern void ODP_Init_aspl_general_s(int nodes, int degree, int* num_degrees, int symmetries);
extern void ODP_Init_aspl_grid(int width, int height, int degree, int* num_degrees);
extern void ODP_Init_aspl_grid_s(int width, int height, int degree, int* num_degrees, int symmetries);
extern void ODP_Set_aspl(void *adjacency, int *diameter, long *sum, double *ASPL);
extern void ODP_Finalize_aspl();

#ifdef MPI_VERSION
#include <mpi.h>
extern void ODP_Init_aspl_mpi_general(int nodes, int degree, int* num_degrees, MPI_Comm comm);
extern void ODP_Init_aspl_mpi_general_s(int nodes, int degree, int* num_degrees, MPI_Comm comm, int symmetries);
extern void ODP_Init_aspl_mpi_grid(int width, int height, int degree, int* num_degrees, MPI_Comm comm);
extern void ODP_Init_aspl_mpi_grid_s(int widht, int height, int degree, int* num_degrees, MPI_Comm comm, int symmetries);
#endif

// APSP utilities
extern void ODP_Set_lbounds_general(int nodes, int degree, int *low_diameter, double *low_ASPL);
extern void ODP_Set_lbounds_grid(int m, int n, int degree, int length, int *low_diameter, double *low_ASPL);
extern void ODP_Set_degrees(int nodes, int lines, int (*edge)[2], int* num_degrees);
extern int  ODP_Get_lines(char *fname);
extern int  ODP_Get_nodes(int lines, int (*edge)[2]);
extern int  ODP_Get_degree(int nodes, int lines, int (*edge)[2]);
extern int  ODP_Get_length(int lines, int height, int (*edge)[2]);
extern bool ODP_Check_general(char *fname);
extern bool ODP_Check_multiple_edges(int lines, int (*edge)[2]);
extern bool ODP_Check_loop(int lines, int (*edge)[2]);
extern void ODP_Read_edge_general(char* fname, int (*edge)[2]);
extern void ODP_Read_edge_grid(char *fname, int *w, int *h, int (*edge)[2]);
extern void ODP_Write_edge_general(int lines, int (*edge)[2], char *fname);
extern void ODP_Write_edge_grid(int lines, int height, int (*edge)[2], char *fname);
extern void ODP_Conv_edge2adjacency_general(int nodes, int lines, int degree, int (*edge)[2], void *adjacency);
extern void ODP_Conv_edge2adjacency_grid(int width, int height, int lines, int degree, int (*edge)[2], void *adjacency);
extern void ODP_Conv_edge2adjacency_general_s(int nodes, int lines, int degree, int (*edge)[2], int symmetries, void *adjacency);
extern void ODP_Conv_edge2adjacency_grid_s(int width, int height, int lines, int degree, int (*edge)[2], int symmetries, int (*adjacency)[degree]);
extern void ODP_Conv_adjacency2edge_general(int nodes, int degree, int *num_degrees, void *adjacency, int (*edge)[2]);
extern void ODP_Conv_adjacency2edge_grid(int width, int height, int degree, int *num_degrees, void *adjacency, int (*edge)[2]);
extern void ODP_Conv_adjacency2edge_general_s(int nodes, int degree,int *num_degrees, void *adjacency, int symmetries, int (*edge)[2]);
extern void ODP_Conv_adjacency2edge_grid_s(int width, int height, int degree,int *num_degrees, void *adjacency, int symmetries, int (*edge)[2]);
extern void ODP_Srand(unsigned int seed);
extern void ODP_Generate_random_general(int nodes, int degree, int (*edge)[2]);
extern void ODP_Generate_random_general_s(int nodes, int degree, int symmetries, int (*edge)[2]);
extern void ODP_Generate_random_grid(int width, int height, int degree, int length, int (*edge)[2]);
extern void ODP_Generate_random_grid_s(int width, int height, int degree, int length, int symmetries, int (*edge)[2]);
extern void ODP_Print_adjacency(int nodes, int degree, int *num_degrees, void *adjacency);
extern void ODP_Print_edge_general(int lines, int (*edge)[2]);
extern void ODP_Print_edge_grid(int lines, int height, int (*edge)[2]);
extern void ODP_Mutate_adjacency_general(int nodes, int degree, int *num_degrees, int (*adjacency)[degree]);
extern void ODP_Mutate_adjacency_general_s(int nodes, int degree, int *num_degrees, int symmetries, int (*adjacency)[degree]);
extern void ODP_Mutate_adjacency_grid(int width, int height, int degree, int *num_degrees, int length, int (*adjacency)[degree]);
extern void ODP_Mutate_adjacency_grid_s(int width, int height, int degree, int *num_degrees, int length, int symmetries, int (*adjacency)[degree]);
extern void ODP_Restore_adjacency_general(int nodes, int degree, int (*adjacency)[degree]);
extern void ODP_Restore_adjacency_general_s(int nodes, int degree, int symmetries, int (*adjacency)[degree]);
extern void ODP_Restore_adjacency_grid(int width, int height, int degree, int (*adjacency)[degree]);
extern void ODP_Restore_adjacency_grid_s(int width, int height, int degree, int symmetries, int (*adjacency)[degree]);
#endif
