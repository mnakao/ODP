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
#define WIDTH(v,h) ((v)/(h))
#define HEIGHT(v,h) ((v)%(h))
#define MAX(x,y) ((x)>(y)?(x):(y))

bool accept(const int nodes, const int current_diameter, const int diameter,
	    const double current_ASPL, const double ASPL, const double temp, const bool enable_ASPL_priority);
bool accept_s(const int nodes, const int current_diameter, const int diameter,
	      const double current_ASPL, const double ASPL, const double temp, const bool enable_ASPL_priority, const int symmetgries);
double get_time();
#endif
