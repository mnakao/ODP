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
#if defined(__ARM_NEON)
#include <arm_neon.h>
#elif !defined(__FUJITSU)
#include <immintrin.h>
#endif
#ifdef _OPENMP
  #include <omp.h>
#endif
#ifdef __NVCC__
  #include <cuda.h>
#endif

typedef struct {
  int u[2], v[2], u_d[2], v_d[2];
} ODP_Restore;

#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define WIDTH(v,h) ((v)/(h))
#define HEIGHT(v,h) ((v)%(h))

#define UINT64_BITS       64
#define ASPL_MATRIX        1
#define ASPL_MATRIX_SAVING 2
#define ASPL_BFS           3
#define NOT_VISITED        0
#define VISITED            1
#define NOT_DEFINED       -1

#if defined(__ARM_NEON) || defined(__FUJITSU)
#define POPCNT(a) __builtin_popcountl(a)
#elif defined(__NVCC__)
#define POPCNT(a) __popcll(a)
#else
#define POPCNT(a) _mm_popcnt_u64(a)
#endif

extern void ODP_Set_lbounds_general(const int nodes, const int degree, int *low_diameter, double *low_ASPL);
extern void ODP_Set_lbounds_grid(const int m, const int n, const int degree, const int length, int *low_diameter, double *low_ASPL);
extern void ODP_Set_degrees(const int nodes, const int lines, int (*edge)[2], int* num_degrees);
extern void ODP_Read_edge_general(const char* fname, int (*edge)[2]);
extern void ODP_Read_edge_grid(const char *fname, int *w, int *h, int (*edge)[2]);
extern void ODP_Conv_edge2adjacency_general(const int nodes, const int lines, const int degree, const int (*edge)[2], int *adjacency);
extern void ODP_Conv_edge2adjacency_grid(const int width, const int height, const int lines, const int degree, const int (*edge)[2], int *adjacency);
extern void ODP_Conv_adjacency2edge_general(const int nodes, const int degree, const int *num_degrees, const int *adjacency, int (*edge)[2]);
extern void ODP_Conv_adjacency2edge_grid(const int width, const int height, const int degree, const int *num_degrees, const int *adjacency, int (*edge)[2]);
extern int ODP_Get_lines(const char* fname);
extern int ODP_Get_nodes(const int lines, const int (*edge)[2]);
extern int ODP_Get_degree(const int nodes, const int lines, const int (*edge)[2]);
extern int ODP_Get_length(const int lines, const int height, const int (*edge)[2]);
#ifdef _OPENMP
extern void ODP_declare_local_frontier(const int nodes);
extern void ODP_free_local_frontier();
#endif
#endif
