#ifndef SA_INCLUDED
#define SA_INCLUDED
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include "odp.h"

#define MAX_FILENAME_LENGTH 256
#define NOT_DEFINED -1
#define ERROR(...) do{fprintf(stderr, __VA_ARGS__); exit(0);}while(0)
#define PRINT_R0(...) do{if(rank == 0){ printf(__VA_ARGS__); }}while(0)
#define MAX(x,y) ((x)>(y)?(x):(y))

bool accept(const int nodes, const int current_diameter, const int diameter,
	    const double current_ASPL, const double ASPL, const double temp, const bool hill_climbing, const bool ASPL_priority);
bool accept_s(const int nodes, const int current_diameter, const int diameter,
	      const double current_ASPL, const double ASPL, const double temp, const bool hill_climbing, const bool ASPL_priority, const int symmetgries);
double get_time();
int get_random(const int max);
int get_degree_index(const int u, const int v, const int u_d, const int nodes,
                     const int degree, const int (*adjacency)[degree]);
bool check_multiple_edges(const int u, const int u_d, const int v,
                          const int degree, const int (*adjacency)[degree]);
void backup_adjacency(const int num, const int u[2], const int u_d[2], const int v[2], const int v_d[2]);
void undo_adjacency(const int num, const int degree, int (*adjacency)[degree]);
#endif
