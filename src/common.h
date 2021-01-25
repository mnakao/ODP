#ifndef COMMON_INCLUDED
#define COMMON_INCLUDED

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include "parameter.h"
#include <immintrin.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#ifdef __NVCC__
  #include <cuda.h>
#endif

#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define UINT64_BITS 64
#define ASPL_NORMAL 1
#define ASPL_SAVING 2
#define NOT_VISITED 0
#define VISITED     1
#define MUTATE_1OPT 1
#define	MUTATE_2OPT 2

#ifdef __NVCC__
#define POPCNT(a) __popcll(a)
#else
#define POPCNT(a) _mm_popcnt_u64(a)
#endif

extern void ODP_Set_lbounds_general(const int nodes, const int degree, int *low_diameter, double *low_ASPL);
extern void ODP_Set_lbounds_grid(const int m, const int n, const int degree, const int length, int *low_diameter, double *low_ASPL);
extern void ODP_Set_degrees(const int nodes, const int lines, int (*edge)[2], int* num_degrees);
extern void ODP_Read_edge_general(const char* fname, int (*edge)[2]);
extern void ODP_Read_edge_grid(const char *fname, int *w, int *h, int (*edge)[2]);
extern void ODP_Conv_edge2adjacency(const int nodes, const int lines, const int degree, const int (*edge)[2], int *adjacency);
extern void ODP_Conv_adjacency2edge(const int nodes, const int degree, const int *num_degrees, const int *adjacency, int (*edge)[2]);
extern int ODP_Get_lines(const char* fname);
extern int ODP_Get_nodes(const int lines, const int (*edge)[2]);
extern int ODP_Get_degree(const int nodes, const int lines, const int (*edge)[2]);
extern int ODP_Get_length(const int lines, const int (*edge)[2], const int height);
#endif
