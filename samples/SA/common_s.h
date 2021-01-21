#ifndef SA_INCLUDED
#define SA_INCLUDED
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include "odp.h"
#define MUTATE_1OPT 1
#define MUTATE_2OPT 2
#define MAX_FILENAME_LENGTH 256
#define ERROR(...) do{fprintf(stderr, __VA_ARGS__); exit(0);}while(0)

void mutate_adjacency_general_s(const int nodes, const int degree, int symmetries, int adjacency[nodes][degree]);
void restore_adjacency(const int degree, int *adjacency);
bool accept_s(const int nodes, const int current_diameter, const int diameter,
	      const double current_ASPL, const double ASPL, const double temp, const int symmetries);
double get_time();
#endif
