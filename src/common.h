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
#include "parameter.h"
#ifdef __AVX2__
#include <immintrin.h>
#endif
#ifdef _OPENMP
  #include <omp.h>
#endif
#ifdef __NVCC__
  #include <cuda.h>
#endif

#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define WIDTH(v, height) (v/height)
#define HEIGHT(v, height) (v%height)
#define UINT64_BITS 64
#define APSP_NORMAL 1
#define APSP_SAVING 2

#ifdef __NVCC__
#define POPCNT(a) __popcll(a)
#else
#define POPCNT(a) _mm_popcnt_u64(a)
#endif

extern void apsp_set_adjacency(const int nodes, const int degree, const int lines,
			      const  int (*edge)[2], int *adjacency);
extern void apsp_set_lbounds_general(const int nodes, const int degree, int *low_diameter,
				     double *low_ASPL);
extern void apsp_set_degrees(const int nodes, const int lines, int (*edge)[2],
			     int* num_degrees);
extern void apsp_read_edge_general(const char* fname, int (*edge)[2]);
extern void apsp_read_edge_grid(const char *fname, int *w, int *h, int (*edge)[2]);
extern void apsp_set_lbounds_grid(const int m, const int n, const int degree, const int length,
				  int *low_diameter, double *low_ASPL);
extern int apsp_get_lines(const char* fname);
extern int apsp_get_nodes(const int lines, int (*edge)[2]);
extern int apsp_get_degree(const int nodes, const int lines, int (*edge)[2]);
extern int apsp_get_length(const int lines, int (*edge)[2], const int height);
#endif
