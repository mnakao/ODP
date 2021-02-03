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

bool accept(const int nodes, const int current_diameter, const int diameter,
	    const double current_ASPL, const double ASPL, const double temp, const bool enable_ASPL_priority);
bool accept_s(const int nodes, const int current_diameter, const int diameter,
	      const double current_ASPL, const double ASPL, const double temp, const bool enable_ASPL_priority, const int symmetgries);
double get_time();
#endif
