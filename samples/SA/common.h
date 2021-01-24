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
#define NOT_VISITED 0
#define VISITED     1
#define MUTATE_1OPT 1
#define MUTATE_2OPT 2
#define ERROR(...) do{fprintf(stderr, __VA_ARGS__); exit(0);}while(0)

void mutate_adjacency_general(const int nodes, const int degree, const int *restrict num_degrees, int adjacency[nodes][degree]);
void mutate_adjacency_general_s(const int nodes, const int degree, const int *restrict num_degrees, const int symmetgries, int adjacency[nodes][degree]);
void mutate_adjacency_grid(const int width, const int height, const int degree, const int *restrict num_degrees, const int length, int (*adjacency)[degree]);
void restore_adjacency(int *adjacency);
bool accept(const int nodes, const int current_diameter, const int diameter,
	    const double current_ASPL, const double ASPL, const double temp);
bool accept_s(const int nodes, const int current_diameter, const int diameter,
	      const double current_ASPL, const double ASPL, const double temp, const int symmetgries);
double get_time();
#endif
