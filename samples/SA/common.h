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
#define ERROR(...) do{fprintf(stderr, __VA_ARGS__); exit(0);}while(0)
#define WIDTH(v, height)  (v/height)
#define HEIGHT(v, height) (v%height)

void mutate_adjacency_general(const int nodes, const int degree, int adjacency[nodes][degree]);
void mutate_adjacency_grid(const int width, const int height, const int degree,
			   const int length, int (*adjacency)[degree]);
void restore_adjacency(const int degree, int *adjacency);
bool accept(const int nodes, const int current_diameter, const int diameter,
	    const double current_ASPL, const double ASPL, const double temp);
double get_time();
#endif
