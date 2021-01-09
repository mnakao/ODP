## Overview
This software is to solve APSP (All Pairs Shortest Paths) for a noweight undirected graph.
You can create following libraries.
* libapsp.a : Serial version
* libapsp_threads.a : Threads version
* libapsp_mpi.a : MPI (Message passing Interface) version
* libapsp_mpi_threads.a: MPI + Threads version
* libapsp_cuda.a: CUDA (Compute Unified Device Architecture) version
* libapsp_mpi_cuda.a: MPI + CUDA version

## Algotherm
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
$ git clone https://github.com/mnakao/APSP.git
$ cd APSP/src
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
$ cd APSP/src
$ make [serial|threads|mpi|mpi_threads|cuda|mpi_cuda|all]
```
 
## How to run sample programs
```
$ cd APSP/sample
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
### APSP=[NORMAL|SAVING]

This library provides two types of APSP algorithms.
One is `NORMAL`, the other is `SAVING`. `SAVING` is a memory-saving version.
By default, `NORMAL` is automatically selected if the amount of memory used is
lower than the value of `MEM_THRESHOLD` in `src/parameter.h`.
This environment variable can specify the algorithm to use.

### APSP_PROFILE=1

Output the performance profile for apsp\*_run().
```
$ APSP=SAVING APSP_PROFILE=1 ./general.x ./graphs/general/n16d4.edges
------ Profile for APSP_RUN ------
Hostname        = kiwi
Initialize Date = Mon Jan  4 23:14:03 2021
Finalize Date   = Mon Jan  4 23:14:03 2021
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
* Hostname : Name of the machine on which the program was run.
* Initialize Date : Time when the apsp\*_init() was executed.
* Finalize Date   : Time when the apsp\*_finalize() was executed.
* Number of Times : Number of times apsp\*_run() was executed.
* Total Time : Total execution time of apsp\*_run()
* Average Time : Average execution time of apsp\*_run()
* Algorithm : NORMAL or SAVING. The parentheses are the types of libraries used. That is, SERIAL, THREADS, MPI, MPI+THREADS, CUDA, or MPI+CUDA.
* Symmetries : When using apsp\*_init(), the value is 1. When using apsp\*_init_s(), the value is symmetries.
* Memory Usage : Amount of memory used in the library.
* Num of Procs : Number of processes used in the library.
* Num of Threads : Number of threads used in the library.

## Basic Function
Note that there are no special functions for the threaded versions.
Thread parallelization is performed automatically depending on the library to be linked.

### Initialize
Perform the initialization process before running APSP.
```
void apsp_init         (int nodes, int degree, int num_degrees[nodes])
void apsp_mpi_init     (int nodes, int degree, int num_degrees[nodes], MPI_Comm comm)
void apsp_cuda_init    (int nodes, int degree, int num_degrees[nodes])
void apsp_mpi_cuda_init(int nodes, int degree, int num_degrees[nodes], MPI_Comm comm)
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree: Degree in a graph.
* [IN] num_degrees : Specify NULL for a regular graph. Or specify the degrees in each vertex for a non-regular graph.
* [IN] comm : MPI communicator.

### Finalize
Release the resources allocated in apsp\*_init().
```
void apsp_finalize()
void apsp_mpi_finalize()
void apsp_cuda_finalize()
void apsp_mpi_cuda_finalize()
```

### Run APSP
Calculate APSP. Note that they must be executed between apsp\*_init() and apsp\*_finalize().
```
void apsp_run         (int adjacency[nodes][degree], int *diameter, long *sum, double *ASPL)
void apsp_mpi_run     (int adjacency[nodes][degree], int *diameter, long *sum, double *ASPL)
void apsp_cuda_run    (int adjacency[nodes][degree], int *diameter, long *sum, double *ASPL)
void apsp_mpi_cuda_run(int adjacency[nodes][degree], int *diameter, long *sum, double *ASPL)
```
* [IN] adjacency : Adjacency matrix of a graph.
* [OUT] diameter : Diameter of a graph.
* [OUT] sum : Total value of the distances between each vertex in a graph.
* [OUT] ASPL : Average shortest path length of a graph (sum = ASPL*(nodes*(nodes-1)/2)).

## Utility
### Read an edge from a file
```
void apsp_read_edge_general(char* fname, int edge[lines][2])
void apsp_read_edge_grid   (char *fname, int *width, int *height, int edge[lines][2])
```
* [IN] fname : File name of a graph.
* [OUT] edge : Edge list of a graph.
* [OUT] width : Width of a grid graph.
* [OUT] height : Height of a grid graph.

### Write an edge to a file
```
void apsp_write_edge_general(int lines, int edge[lines][2], char *fname)
void apsp_write_edge_grid   (int lines, int height, int edge[lines][2], char *fname)
```
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.
* [IN] height : Height of a grid graph.
* [OUT] fname : File name of a graph.

### Convert an edge list to an adjacency matrix
```
void apsp_conv_edge2adjacency(int nodes, int lines, int edge[lines][2], int adjacency[nodes][degree])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.
* [OUT] adjacency : Adjacency matrix of a graph.

### Convert an adjacency matrix to an edge list
```
void apsp_conv_adjacency2edge(int nodes, int degree, int num_degrees[nodes], int adjacency[nodes][degree], int edge[lines][2])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree : Degree in a graph.
* [IN] num_degrees : Specify NULL for a regular graph. If not, specify the degrees for each vertex.
* [IN] adjacency : Adjacency matrix of a graph.
* [OUT] edge : Edge list of a graph.

### Set theoretical lower bounds
```
void apsp_set_lbounds_general(int nodes, int degree, int *low_diameter, double *low_ASPL)
void apsp_set_lbounds_grid   (int width, int height, int degree, int length, int *low_diameter, double *low_ASPL)
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
void apsp_set_degrees(int nodes, int lines, int edge[lines][2], int num_degrees[nodes])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.
* [OUT] num_degrees : Degree in each vertex.

### Set the seed value for random number generation
This function is executed before apsp_generate_random\*() and apsp_mutate_adjencency\*() that uses random numbers internally.
```
void apsp_srand(unsigned int seed)
```
* [IN] seed : Seed for random.

### Generate a random graph
Generate a regular graph with randomly connected vertices. Note that the graph may contain multiple edges and loops.
```
void apsp_generate_random_general(int nodes, int degree, int edge[lines][2])
void apsp_generate_random_grid   (int width, int height, int degree, int length, int edge[lines][2])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree : Degree in a graph.
* [IN] width : Width of a grid graph.
* [IN] height : Height of a grid graph.
* [IN] length : Maximum length of a grid graph.
* [OUT] edge : Edge list of a graph.

### Get the number of lines in a file
```
int apsp_get_lines(char* fname)
```
* [RETURN] : Nnumber of lines in a file.
* [IN] fname : File name of a graph.

### Get the number of nodes in a graph
```
int apsp_get_nodes(int lines, int edge[lines][2])
```
* [RETURN] : Nnumber of nodes in an edge list.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.

### Get a degree in a graph
```
int apsp_get_degree(int nodes, int lines, int edge[lines][2])
```
* [RETURN] : Degree in an edge list.
* [IN] nodes : Number of nodes in a graph.
* [IN] edge : Edge list of a graph.

### Get a maximum length for a grid graph
```
int apsp_get_length(int lines, int edge[lines][2], int height)
```
* [RETURN] : Maximum length in an edge list.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.
* [IN] height : Height of a grid graph.

### Check if an input file is a general graph
```
bool apsp_check_general(char *fname)
```
* [RETUREN] : When an input is a general graph, it returns true.
* [IN] fname : File name of a graph.

### Check if a graph has multiple edges
```
bool apsp_check_multiple_edges(int lines, int edge[lines][2])
```
* [RETURN] : If a graph has multiple edges, it returns true.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.

### Check if a graph has a self-loop
```
bool apsp_check_loop(int lines, int edge[lines][2])
```
* [RETURN] : If a graph has a self-loop, it returns true.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.

### Shortcut functions
Return other values just by specifying fname and comm.
Please see `samples/simple_general.c`.

```
void apsp_all_run_general         (char *fname, int *nodes, int *degree, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_run_grid            (char *fname, int *width, int *height, int *degree, int *length, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_mpi_run_general     (char *fname, MPI_Comm comm, int *nodes, int *degree, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_mpi_run_grid        (char *fname, MPI_Comm comm, int *width, int *height, int *degree, int *length, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_cuda_run_general    (char *fname, int *nodes, int *degree, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_cuda_run_grid       (char *fname, int *width, int *height, int *degree, int *length, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_mpi_cuda_run_general(char *fname, MPI_Comm comm, int *nodes, int *degree, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_mpi_cuda_run_grid   (char *fname, MPI_Comm comm, int *width, int *height, int *degree, int *length, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
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

## Function for a graph with symmetry
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

### Initialize
These functions can be used instead of the apsp\*_init() for only a general graph with symmetry.

```
void apsp_init_s         (int nodes, int degree, int num_degrees[nodes], int symmetries)
void apsp_mpi_init_s     (int nodes, int degree, int num_degrees[nodes], MPI_Comm comm, int symmetries)
void apsp_cuda_init_s    (int nodes, int degree, int num_degrees[nodes], int symmetries)
void apsp_mpi_cuda_init_s(int nodes, int degree, int num_degrees[nodes], MPI_Comm comm, int symmetries)
```
* [IN] nodes : Number of nodes in a graph.
* [IN] degree: Degree in a graph.
* [IN] num_degrees : Specify NULL for a regular graph. If not, specify the degrees for each vertex.
* [IN] comm : MPI communicator.
* [IN] symmetries : Numer of symmetries in a graph. This value must be a divisor of nodes. If it is 1, it works the same as apsp\*_init().

Note that the apsp\*_finalize() can be used in common.

### Convert an edge list to an adjacency matrix
```
void apsp_conv_edge2adjacency_s(int nodes, int lines, int edge[lines][2], int symmetries, int adjacency[nodes/symmetries][degree])
```
* [IN] nodes : Number of nodes in a graph.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list of a graph.
* [IN] symmetries : Numer of symmetries in a graph. This value must be a divisor of nodes. If it is 1, it works the same as apsp_conv_edge2adjacency().
* [OUT] adjacency : Adjacency matrix of a graph.

## Note
The software also supports non-regular graphs, but usage of memory may be not efficient.
Because the format of the adjacency matrix is `int adjacency[nodes][degree]`, which is commanly used in regular and non-regular graphs.
This wastes memory when the number of edges that each vertex has varies greatly.

## Acknowledgment
This software is inherited from https://github.com/ryuhei-mori/graph_ASPL.git.
