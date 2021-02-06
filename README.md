## Overview
This software is to solve Order/Degree Problem for a noweight undirected graph.
You can create following libraries.
* libodp.a : Serial version
* libodp_threads.a : Threads version
* libodp_mpi.a : MPI (Message passing Interface) version
* libodp_mpi_threads.a: MPI + Threads version
* libodp_cuda.a: CUDA (Compute Unified Device Architecture) version
* libodp_mpi_cuda.a: MPI + CUDA version

## Algotherm to obtain ASPL (Average Shortest Path Length)
Please see the paper (Open Access).
* https://dl.acm.org/doi/10.1145/3368474.3368478

@inproceedings{10.1145/3368474.3368478,
  author = {Nakao, Masahiro and Murai, Hitoshi and Sato, Mitsuhisa},
  title = {Parallelization of All-Pairs-Shortest-Path Algorithms in Unweighted Graph},
  year = {2020},
  isbn = {9781450372367},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  url = {https://doi.org/10.1145/3368474.3368478},
  doi = {10.1145/3368474.3368478},
  booktitle = {Proceedings of the International Conference on High Performance Computing in Asia-Pacific Region},
  pages = {63–72},
  numpages = {10},
  keywords = {network, hybrid parallelization, graph theory, GPU},
  location = {Fukuoka, Japan},
  series = {HPCAsia2020}
}

_When you write a paper using this software, please refer to the paper._

## Quick start
```
$ git clone https://github.com/mnakao/ODP.git
$ cd ODP/src
$ make
$ cd ../samples
$ make
$ ./general.x ./graphs/general/n16d4.edges
Nodes = 16, Degrees = 4
Diameter     = 3
Diameter Gap = 1 (3 - 2)
ASPL         = 1.9166666667 (230/120)
ASPL Gap     = 0.1833333333 (1.9166666667 - 1.7333333333)
```

## How to create libraries
```
$ cd ODP/src
$ make [serial|threads|mpi|mpi_threads|cuda|mpi_cuda|all]
```
 
## How to run sample programs
```
$ cd ODP/sample
$ make [serial|threads|mpi|mpi_threads|cuda|mpi_cuda|all]
$ ./general.x ./graphs/general/n16d4.edges
Nodes = 16, Degrees = 4
Diameter     = 3
Diameter Gap = 1 (3 - 2)
ASPL         = 1.9166666667 (230/120)
ASPL Gap     = 0.1833333333 (1.9166666667 - 1.7333333333)
```

## Graph file format
* A graph file must be in an edge list format compatible with the NetworkX's read_edgelist function. Edge attributes will be ignored if exist.
* For a general graph, each vertex name must be an integer starting from zero.
* For a grid graph, each vertex name must be a comma-separated string like "x,y" (no quotes, no spaces). x and y are integers starting from zero, which represent the coordinate of the vertex.
* Please see sample graphs in `./samples/graphs/` or http://research.nii.ac.jp/graphgolf/submit.html

## Environment variable
### ODP_ASPL=[NORMAL|SAVING]

This library provides two types of algorithms to obtain ASPL.
One is `NORMAL`, the other is `SAVING`. `SAVING` is a memory-saving version of `NORMAL`.
By default, `NORMAL` is automatically selected if the amount of memory used is
lower than the value of `MEM_THRESHOLD` in `src/parameter.h`.
This environment variable can specify the algorithm to use.

### ODP_PROFILE=1

Output the performance profile for ODP_Set_aspl(). 
This profile is output when ODP_Finalize_aspl() is executed.
```
$ ODP_ASPL=SAVING ODP_PROFILE=1 ./general.x ./graphs/general/n16d4.edges
------ Profile for SET_ASPL ------
Date            = Mon Jan  4 23:14:03 2021
Hostname        = kiwi
Number of Times = 1
Total Time      = 0.000005 sec.
Average Time    = 0.000005 sec.
Algorithm       = SAVING (SERIAL)
Symmetries      = 1
Memory Usage    = 0.006 MB
Num of Procs    = 1
Num of Threads  = 1
--------- End of Profile ---------
Nodes = 16, Degrees = 4
Diameter     = 3
Diameter Gap = 1 (3 - 2)
ASPL         = 1.9166666667 (230/120)
ASPL Gap     = 0.1833333333 (1.9166666667 - 1.7333333333)
```

The meaning of each item in the profile is as follows.
* Date : The time when the profile was output.
* Hostname : Name of the machine on which the program was run.
* Number of Times : Number of times ODP_Set_aspl() was executed.
* Total Time : Total execution time of ODP_Set_aspl()
* Average Time : Average execution time of ODP_Set_aspl()
* Algorithm : NORMAL or SAVING. The parentheses are the types of libraries used. That is, SERIAL, THREADS, MPI, MPI+THREADS, CUDA, or MPI+CUDA.
* Symmetries : When using ODP_Init_aspl\*_s(), the value is `symmetries`. Otherwise, the value is 1.
* Memory Usage : Amount of memory used in the library.
* Num of Procs : Number of processes used in the library.
* Num of Threads : Number of threads used in the library.

## Basic Function
Note that there are no special functions for the threaded versions.
Thread parallelization is performed automatically depending on the library to be linked.

### Initialize
Perform the initialization process before executing ODP_Set_aspl().
```
void ODP_Init_aspl_general         (int nodes, int degree, int num_degrees[nodes])
void ODP_Init_aspl_grid            (int width, int height, int degree, int num_degrees[nodes])
void ODP_Init_aspl_mpi_general     (int nodes, int degree, int num_degrees[nodes], MPI_Comm comm)
void ODP_Init_aspl_mpi_grid        (int width, int height, int degree, int num_degrees[nodes], MPI_Comm comm)
void ODP_Init_aspl_cuda_general    (int nodes, int degree, int num_degrees[nodes])
void ODP_Init_aspl_cuda_grid       (int width, int height, int degree, int num_degrees[nodes])
void ODP_Init_aspl_mpi_cuda_general(int nodes, int degree, int num_degrees[nodes], MPI_Comm comm)
void ODP_Init_aspl_mpi_cuda_grid   (int width, int height, int degree, int num_degrees[nodes], MPI_Comm comm)
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree: Degree in a graph.
* [IN] num_degrees : Specify NULL for a regular graph. Or specify the degrees in each vertex for a non-regular graph.
* [IN] width : Width of a grid graph.
* [IN] height : Height of a grid graph.
* [IN] comm : MPI communicator.

### Set diameter, sum, and ASPL
Set diameter, sum, and ASPL. Note that these functions must be executed between ODP_Init_aspl\*() and ODP_Finalize_aspl().
In the case of an unconnected graph, INT_MAX, LONG_MAX, and DBL_MAX are assigned to the values of diameter, sum, and ASPL, respectively.
```
void ODP_Set_aspl(int adjacency[nodes][degree], int *diameter, long *sum, double *ASPL)
```
* [IN] adjacency : Adjacency matrix of a graph.
* [OUT] diameter : Diameter of a graph.
* [OUT] sum : Total value of the distances between each vertex in a graph.
* [OUT] ASPL : Average shortest path length of a graph (sum = ASPL*(nodes*(nodes-1)/2)).

### Finalize
Release the resources allocated in ODP_Init_aspl\*().
```
void ODP_Finalize_aspl()
```

## Utility
### Read an edge from a file
```
void ODP_Read_edge_general(char* fname, int edge[lines][2])
void ODP_Read_edge_grid   (char *fname, int *width, int *height, int edge[lines][2])
```
* [IN] fname : File name of a graph.
* [OUT] edge : Edge list of a graph.
* [OUT] width : Width of a grid graph.
* [OUT] height : Height of a grid graph.

### Write an edge to a file
```
void ODP_Write_edge_general(int lines, int edge[lines][2], char *fname)
void ODP_Write_edge_grid   (int lines, int height, int edge[lines][2], char *fname)
```
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.
* [IN] height : Height of a grid graph.
* [OUT] fname : File name of a graph.

### Print an adjacency matrix
```
void ODP_Print_adjacency(int nodes, int degree, int num_degrees[nodes], int adjacency[nodes][degree])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree : Degree in a graph.
* [IN] num_degrees : Specify NULL for a regular graph. If not, specify the degrees for each vertex.
* [IN] adjacency : Adjacency matrix of a graph.

### Print an edge list
```
void ODP_Print_edge_general(int lines, int edge[lines][2])
void ODP_Print_edge_grid(int lines, int height, int edge[lines][2])
```
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.
* [IN] height : Height of a grid graph.

### Convert an edge list to an adjacency matrix
```
void ODP_Conv_edge2adjacency(int nodes, int lines, int degree, int edge[lines][2], int adjacency[nodes][degree])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] lines : Number of lines in an edge list.
* [IN] degree : Degree in a graph.
* [IN] edge : Edge list of a graph.
* [OUT] adjacency : Adjacency matrix of a graph.

### Convert an adjacency matrix to an edge list
```
void ODP_Conv_adjacency2edge(int nodes, int degree, int num_degrees[nodes], int adjacency[nodes][degree], int edge[lines][2])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree : Degree in a graph.
* [IN] num_degrees : Specify NULL for a regular graph. If not, specify the degrees for each vertex.
* [IN] adjacency : Adjacency matrix of a graph.
* [OUT] edge : Edge list of a graph.

### Set theoretical lower bounds
```
void ODP_Set_lbounds_general(int nodes, int degree, int *low_diameter, double *low_ASPL)
void ODP_Set_lbounds_grid   (int width, int height, int degree, int length, int *low_diameter, double *low_ASPL)
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree : Degree in a graph.
* [OUT] low_diameter : Theoretical lower bound of diameter in a graph.
* [OUT] low_ASPL : Theoretical lower bound of ASPL in a graph.
* [IN] width : Width of a grid graph.
* [IN] height : Height of a grid graph.
* [IN] length : Maximum length of a grid graph.

### Set degrees for a non-regular graph
```
void ODP_Set_degrees(int nodes, int lines, int edge[lines][2], int num_degrees[nodes])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.
* [OUT] num_degrees : Degree in each vertex.

### Seed for a random number
Since a random number is used in ODP_Generate_random_\*() and ODP_Mutate_adjacency_\*(), ODP_Srand() must be executed before those functions.
```
void ODP_Srand(unsigned int seed)
```
* [IN] seed : Seed for random.

### Generate a random graph
Generate a regular graph with randomly connected vertices. Note that the graph may contain multiple edges and loops.
```
void ODP_Generate_random_general(int nodes, int degree, int edge[lines][2])
void ODP_Generate_random_grid   (int width, int height, int degree, int length, int edge[lines][2])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree : Degree in a graph.
* [IN] width : Width of a grid graph.
* [IN] height : Height of a grid graph.
* [IN] length : Maximum length of a grid graph.
* [OUT] edge : Edge list of a graph.

### Mutate an adjacency matrix
Mutate an adjacency matrix slightly. Specifically, the operation equivalent to the 2-opt method is performed. 
```
void ODP_Mutate_adjacency_general(int nodes, int degree, int num_degrees[nodes], int adjacency[nodes][degree])
void ODP_Mutate_adjacency_grid(int width, int height, int degree, int num_degrees[nodes], int length, int adjacency[nodes][degree])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree : Degree in a graph.
* [IN] num_degrees : Degree in each vertex.
* [IN] width : Width of a grid graph.
* [IN] height : Height of a grid graph.
* [IN] length : Maximum length of a grid graph.
* [OUT] adjacency : Adjacency matrix of a graph.

### Restore an adjacency matrix
Undo before being modified by ODP_Mutate_adjacency_\*().
Only if there is no change in the adjacency matrix between ODP_Mutate_adjacency_\*() and ODP_Restore_adjacency_\*(). 
```
void ODP_Restore_adjacency_general(int nodes, int degree, int adjacency[nodes][degree])
void ODP_Restore_adjacency_grid(int width, int height, int degree, int adjacency[nodes][degree])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree : Degree in a graph.
* [IN] width : Width of a grid graph.
* [IN] height : Height of a grid graph.
* [OUT] adjacency : Adjacency matrix of a graph.

### Get the number of lines in a file
```
int ODP_Get_lines(char* fname)
```
* [RETURN] : Number of lines in a file.
* [IN] fname : File name of a graph.

### Get the number of nodes in a graph
```
int ODP_Get_nodes(int lines, int edge[lines][2])
```
* [RETURN] : Number of nodes in an edge list.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.

### Get a degree in a graph
```
int ODP_Get_degree(int nodes, int lines, int edge[lines][2])
```
* [RETURN] : Degree in an edge list.
* [IN] nodes : Number of nodes in a graph.
* [IN] edge : Edge list of a graph.

### Get a maximum length for a grid graph
```
int ODP_Get_length(int lines, int edge[lines][2], int height)
```
* [RETURN] : Maximum length in an edge list.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.
* [IN] height : Height of a grid graph.

### Check if an input file is a general graph
```
bool ODP_Check_general(char *fname)
```
* [RETURN] : When an input is a general graph, it returns true.
* [IN] fname : File name of a graph.

### Check if a graph has multiple edges
```
bool ODP_Check_multiple_edges(int lines, int edge[lines][2])
```
* [RETURN] : If a graph has multiple edges, it returns true.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.

### Check if a graph has a self-loop
```
bool ODP_Check_loop(int lines, int edge[lines][2])
```
* [RETURN] : If a graph has a self-loop, it returns true.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.

### Shortcut functions
Return other values just by specifying fname and comm.
Please see `samples/simple_general.c`.

```
void ODP_Set_aspl_general         (char *fname, int *nodes, int *degree, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void ODP_Set_aspl_grid            (char *fname, int *width, int *height, int *degree, int *length, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void ODP_Set_aspl_mpi_general     (char *fname, MPI_Comm comm, int *nodes, int *degree, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void ODP_Set_aspl_mpi_grid        (char *fname, MPI_Comm comm, int *width, int *height, int *degree, int *length, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void ODP_Set_aspl_cuda_general    (char *fname, int *nodes, int *degree, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void ODP_Set_aspl_cuda_grid       (char *fname, int *width, int *height, int *degree, int *length, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void ODP_Set_aspl_mpi_cuda_general(char *fname, MPI_Comm comm, int *nodes, int *degree, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void ODP_Set_aspl_mpi_cuda_grid   (char *fname, MPI_Comm comm, int *width, int *height, int *degree, int *length, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
```
* [IN] fname : File name of a graph.
* [IN] comm : MPI communicator.
* [OUT] nodes : Number of nodes in a graph.
* [OUT] degree : Degree in a graph.
* [OUT] low_diameter : Theoretical lower bound of diameter in a graph.
* [OUT] low_ASPL : Theoretical lower bound of ASPL in a graph.
* [OUT] diameter : Diameter of a graph.
* [OUT] sum : Total value of the distances between each vertex in a graph.
* [OUT] ASPL : Average shortest path length of a graph (sum = ASPL*(nodes*(nodes-1)/2)).
* [OUT] width : Width of a grid graph.
* [OUT] height : Height of a grid graph.
* [OUT] length : Maximum length of a grid graph.

## Function for a general graph with symmetry
Symmetry in this software means that when the vertices are arranged on a circle,
the graph when rotated by `360/symmetries` degrees and the graph before rotation match.

![](./misc/img/symmetry.png)

The above image is an example of a graph with (nodes, degree, symmetries) = (24, 4, 4).
This image is created from `samples/graphs/general/n24d4g4.png` and `samples/graphs/general/n24d4g4.edges`.
It also shows the adjacency matrix created from the edge list.
The adjacency matrix can be divided into four groups (`= symmetries`).
The values on the 1st row are 19, 9, and 2.
It means that the vertex number 0 has three edges, 0-19, 0-9, and 0-2.
The edges plus 6 (`= nodes/symmetries`) matches the 1st row in the next group (in line 7).
2 + 6 = 8 and 9 + 6 = 15.
Here, 19 + 6 = 25, but the number of nodes is 24, so it goes around and becomes 25 - 24 = 1.
This rule holds for all groups.

The remaining elements can be calculated from the first (nodes/symmetries) lines of the adjacency matrix.
Thus, the size of the adjacency matrix is `int adjacency[nodes/symmetries][degree]`.
However, since the edge list is used for input and output, the size is `int edge[lines][2]`, which is the same as a normal graph.

### Initialize
These functions can be used instead of the ODP_Init_aspl\*().

```
void ODP_Init_aspl_general_s         (int nodes, int degree, int num_degrees[nodes/symmetries], int symmetries)
void ODP_Init_aspl_grid_s            (int width, int height, int degree, int num_degrees[nodes/symmetries], int symmetries)
void ODP_Init_aspl_mpi_general_s     (int nodes, int degree, int num_degrees[nodes/symmetries], MPI_Comm comm, int symmetries)
void ODP_Init_aspl_mpi_grid_s        (int width, int height, int degree, int num_degrees[nodes/symmetries], MPI_Comm comm, int symmetries)
void ODP_Init_aspl_cuda_general_s    (int nodes, int degree, int num_degrees[nodes/symmetries], int symmetries)
void ODP_Init_aspl_cuda_grid_s       (int width, int height, int degree, int num_degrees[nodes/symmetries], int symmetries)
void ODP_Init_aspl_mpi_cuda_general_s(int nodes, int degree, int num_degrees[nodes/symmetries], MPI_Comm comm, int symmetries)
void ODP_Init_aspl_mpi_cuda_grid_s   (int width, int height, int degree, int num_degrees[nodes/symmetries], MPI_Comm comm, int symmetries)
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree: Degree in a graph.
* [IN] num_degrees : Specify NULL for a regular graph. If not, specify the degrees for each vertex.
* [IN] width : Width of a grid graph.
* [IN] height : Height of a grid graph.
* [IN] comm : MPI communicator.
* [IN] symmetries : Numer of symmetries in a graph. This value must be a divisor of nodes. If it is 1, it works the same as ODP_Init_aspl\*().

Note that the ODP_Set_aspl() and ODP_Finalize_aspl() can be used in common.

### Convert an edge list to an adjacency matrix
```
void ODP_Conv_edge2adjacency_general_s(int nodes, int lines, int degree, int edge[lines][2], int symmetries, int adjacency[nodes/symmetries][degree])
void ODP_Conv_edge2adjacency_grid_s(int width, int height, int lines, int degree, int edge[lines][2], int symmetries, int adjacency[nodes/symmetries][degree])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] lines : Number of lines in an edge list.
* [IN] degree: Degree in a graph.
* [IN] edge : Edge list of a graph.
* [IN] symmetries : Numer of symmetries in a graph. This value must be a divisor of nodes. If it is 1, it works the same as ODP_Conv_edge2adjacency().
* [OUT] adjacency : Adjacency matrix of a graph.

### Convert an adjacency matrix to an edge list
```
void ODP_Conv_adjacency2edge_general_s(int nodes, int degree, int num_degrees[nodes/symmetries], int adjacency[nodes/symmetries][degree], int symmetries, int edge[lines][2])
void ODP_Conv_adjacency2edge_grid_s(int width, int height, int degree, int num_degrees[nodes/symmetries], int adjacency[nodes/symmetries][degree], int symmetries, int edge[lines][2])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree: Degree in a graph.
* [IN] num_degrees : Specify NULL for a regular graph. If not, specify the degrees for each vertex.
* [IN] adjacency : Adjacency matrix of a graph.
* [IN] symmetries : Numer of symmetries in a graph. This value must be a divisor of nodes. If it is 1, it works the same as ODP_Conv_adjacency2edge().
* [IN] width : Width of a grid graph.
* [IN] height : Height of a grid graph.
* [OUT] edge : Edge list of a graph.

### Generate a random graph
```
void ODP_Generate_random_general_s(int nodes, int degree, unsigned int seed, int symmetries, int edge[lines][2])
void ODP_Generate_random_grid_s(int width, int height, int degree, int length, unsigned int seed, int symmetries, int edge[lines][2])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree: Degree in a graph.
* [IN] seed : Seed for random.
* [IN] symmetries : Numer of symmetries in a graph. This value must be a divisor of nodes. If it is 1, it works the same as ODP_Generate_random_general().
* [IN] width : Width of a grid graph.
* [IN] height : Height of a grid graph.
* [IN] length : Maximum length of a grid graph.
* [OUT] edge : Edge list of a graph.

## Note
The software also supports non-regular graphs, but usage of memory may be not efficient.
Because the format of the adjacency matrix is `int adjacency[nodes][degree]`, which is commanly used in regular and non-regular graphs.
This wastes memory when the number of edges that each vertex has varies greatly.

## Acknowledgment
This software is inherited from https://github.com/ryuhei-mori/graph_ASPL.git.
