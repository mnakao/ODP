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

## How to create libraries
```
$ cd src
$ make [serial|threads|mpi|mpi_threads|cuda|mpi_cuda|all]
```
 
## How to run sample programs
```
$ cd sample
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

Output the performance profile.
```
$ APSP=SAVING APSP_PROFILE=1 ./general.x ./graphs/general/n16d4.edges
------ Start of Profile ------
Hostname       = Pele.local
Start Date     = Sun Dec 27 22:17:15 2020
End Date       = Sun Dec 27 22:17:15 2020
Elapsed Time   = 0.000980 sec.
Algorithm      = SAVING (SERIAL)
Symmetry       = 1
Memory Usage   = 0.016 MB
------  End of Profile  ------
Nodes = 16, Degrees = 4
Diameter     = 3
Diameter Gap = 1 (3 - 2)
ASPL         = 1.9166666667 (230/120)
ASPL Gap     = 0.1833333333 (1.9166666667 - 1.7333333333)
```

## Functions
Note that there are no special functions for the threaded versions.
Thread parallelization is performed automatically depending on the library to be linked.

### Initialize
Perform the initialization process before running APSP.
```
void apsp_init         (int nodes, int degree, int num_degrees[nodes])
void apsp_cuda_init    (int nodes, int degree, int num_degrees[nodes])
void apsp_mpi_init     (int nodes, int degree, int num_degrees[nodes], MPI_Comm comm)
void apsp_mpi_cuda_init(int nodes, int degree, int num_degrees[nodes], MPI_Comm comm)
```
* [IN] nodes : Number of nodes in a graph
* [IN] degree: Degree in a graph
* [IN] num_degrees[nodes] : Specify NULL for a regular graph. Or specify the degrees in each vertex for a non-regular graph.
* [IN] comm : MPI communicator

### Finalize
Release the resources allocated in apsp_\*init().
```
void apsp_finalize()
void apsp_cuda_finalize()
void apsp_mpi_finalize()
void apsp_mpi_cuda_finalize()
```

### Run APSP
Calculate APSP. Note that they must be executed between apsp_\*init() and apsp_\*finalize().
```
void apsp_run         (int adjacency[nodes][degree], int *diameter, long *sum, double *ASPL)
void apsp_cuda_run    (int adjacency[nodes][degree], int *diameter, long *sum, double *ASPL)
void apsp_mpi_run     (int adjacency[nodes][degree], int *diameter, long *sum, double *ASPL)
void apsp_mpi_cuda_run(int adjacency[nodes][degree], int *diameter, long *sum, double *ASPL)
```
* [IN] adjacency[nodes][degree] : Adjacency matrix of a graph
* [OUT] diameter : Diameter of a graph
* [OUT] sum : Total value of the distances between each vertex in a graph
* [OUT] ASPL : Average shortest path length of a graph (sum = ASPL*(nodes*(nodes-1)/2))

### Calculate theoretical lower bounds
```
void apsp_set_lbounds_general(int nodes, int degree, int *low_diameter, double *low_ASPL)
void apsp_set_lbounds_grid   (int width, int height, int degree, int length, int *low_diameter, double *low_ASPL)
```
* [IN] nodes : Number of nodes in a graph
* [IN] degree : Degree in a graph
* [OUT] low_diameter : Theoretical lower bound of diameter in a graph
* [OUT] low_ASPL : Theoretical lower bound of ASPL in a graph
* [IN] width : Width of a grid graph
* [IN] height : Height of a grid graph
* [IN] length : Maximum length of a grid graph

### Set edge from a file
```
void apsp_set_edge_general(char* fname, int (*edge)[2])
void apsp_set_edge_grid   (char *fname, int *width, int *height, int (*edge)[2])
```
* [IN] fname : File name of a graph
* [OUT] edge : Edge list of a graph
* [OUT] width : Width of a grid graph
* [OUT] height : Height of a grid graph

### Set adjacency matrix from an edge list
```
void apsp_set_adjacency(int nodes, int degree, int lines, int edge[lines][2], int adjacency[nodes][degree])
```
* [IN] nodes : Number of nodes in a graph
* [IN] degree : Degree in a graph
* [IN] lines : Number of lines in an edge list
* [IN] edge[lines][2] : Edge list of a graph
* [OUT] adjacency[nodes][degree] : Adjacency matrix of a graph

### Set degrees for a non-regular graph
```
void apsp_set_degrees(int nodes, int lines, int edge[lines][2], int num_degrees[nodes])
```
* [IN] nodes : Number of nodes in a graph
* [IN] lines : Number of lines in an edge list
* [IN] edge[lines][2] :	Edge list of a graph
* [OUT] num_degrees[nodes] : Degree in each vertex

### Find out number of lines in a file
```
int apsp_get_lines(char* fname)
```
* [RETURN] : Nnumber of lines in a file
* [IN] fname : File name of a graph

### Find out number of nodes in a graph
```
int apsp_get_nodes(int lines, int edge[lines][2])
```
* [RETURN] : Nnumber of nodes in an edge list
* [IN] lines : Number of lines in an edge list
* [IN] edge[lines][2] : Edge list of a graph

### Find out degree in a graph
```
int apsp_get_degree(int nodes, int lines, int edge[lines][2])
```
* [RETURN] : Degree in an edge list
* [IN] nodes : Number of nodes in a graph
* [IN] edge[lines][2] : Edge list of a graph

### Find out maximum length for a grid graph
```
int apsp_get_length(int lines, int edge[lines][2], int height)
```
* [RETURN] : Maximum length in an edge list
* [IN] lines : Number of lines in an edge list
* [IN] edge[lines][2] :	Edge list of a graph
* [IN] height : Height of a grid graph

### Find out whether input file is a general graph
```
bool apsp_check_general(char *fname)
```
* [RETUREN] : When an input is a general graph, it returns true.
* [IN] fname : File name of a graph

### Find out whether a graph has duplicated edges
```
bool apsp_check_duplicated_edge(int lines, int edge[lines][2])
```
* [RETURN] : If a graph has duplicated edge, it returns true
* [IN] lines : Number of lines in an edge list
* [IN] edge[lines][2] : Edge list of a graph

### Find out whether a graph has a self-loop
```
bool apsp_check_loop(int lines, int edge[lines][2])
```
* [RETURN] : If a graph has a self-loop, it returns true
* [IN] lines : Number of lines in an edge list
* [IN] edge[lines][2] : Edge list of a graph

### Shortcut functions
Return other values just by specifying fname and comm.
Please see `samples/simple_general.c`.

```
void apsp_all_run_general(char *fname, int *nodes, int *degree, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_run_grid   (char *fname, int *width, int *height, int *degree, int *length, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_cuda_run_general(char *fname, int *nodes, int *degree, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_cuda_run_grid(char *fname, int *width, int *height, int *degree, int *length, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_mpi_run_general(char *fname, MPI_Comm comm, int *nodes, int *degree, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_mpi_run_grid(char *fname, MPI_Comm comm, int *width, int *height, int *degree, int *length, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_mpi_cuda_run_general(char *fname, MPI_Comm comm, int *nodes, int *degree, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
void apsp_all_mpi_cuda_run_grid(char *fname, MPI_Comm comm, int *width, int *height, int *degree, int *length, int *low_diameter, double *low_ASPL, int *diameter, long *sum, double *ASPL)
```
* [IN] fname : File name of a graph
* [IN] comm : MPI communicator
* [OUT] nodes : Number of nodes in a graph
* [OUT] degree : Degree in a graph
* [OUT] low_diameter : Theoretical lower bound of diameter in a graph
* [OUT] low_ASPL : Theoretical lower bound of ASPL in a graph
* [OUT] diameter : Diameter of a graph
* [OUT] sum : Total value of the distances between each vertex in a graph
* [OUT] ASPL : Average shortest path length of a graph (sum = ASPL*(nodes*(nodes-1)/2))
* [OUT] width : Width of a grid graph
* [OUT] height : Height of a grid graph
* [OUT] length : Maximum length of a grid graph

### Initialize for a graph with symmetry
These functions can be used instead of the apsp_init\* functions for only a general graph with symmetry.

```
void apsp_init_s         (int nodes, int degree, int num_degrees[nodes], int groups)
void apsp_cuda_init_s    (int nodes, int degree, int num_degrees[nodes], int groups)
void apsp_mpi_init_s     (int nodes, int degree, int num_degrees[nodes], MPI_Comm comm, int groups)
void apsp_mpi_cuda_init_s(int nodes, int degree, int num_degrees[nodes], MPI_Comm comm, int groups)
```
* [IN] nodes : Number of nodes in a graph
* [IN] degree: Degree in a graph
* [IN] num_degrees[nodes] : Specify NULL for a regular graph. If not, specify the degrees for each vertex.
* [IN] comm : MPI communicator
* [IN] groups : Numer of groups in a graph with symmetry. This value must be a divisor of nodes.

Symmetry in this software means that when the vertices are arranged on a circle,
the graph when rotated by `360/groups` degrees and the graph before rotation match.
Please see `samples/graphs/general/n24d4g4.edges` and `samples/graphs/general/n24d4g4.png` as an example of (nodes, degree, groups) = (24, 4, 4).

![](./misc/img/symmetry.png)

This image is a `samples/graphs/general/n24d4g4.edges`, and is divided into four ( = the value of `groups`).
For example, the values on the 1st row are 22 and 15.
The number of them plus 6 ( `= nodes/groups`) matches the 1st row in the next group (in line 10).
Here, 22 + 6 = 28, but the number of nodes is 24, so it goes around and becomes 28 - 24 = 6.

There is also an exception if the value of `groups` are even.
The edge with vertex numbers 4 and 16 in the 3rd line is the diameter across the center of the circle.
Specifically, the value obtained by subtracting the vertex numbers is half the number of nodes.
The difference between the 3rd and 12th lines is exactly 6, but the difference between the 12th and 21st lines is not 6.
This is because if the difference between the 3rd and 21st lines is 12, the lines will overlap.
In such cases, the edge list is devided into two large groups.
One group consists of the 1st to `groups/2`th group, and another group consists of the remaining groups (`groups/2+1`th to the last group).
Within each group, the difference between the vertices should be `nodes/groups`.

## Note
The software also supports non-regular graphs, but usage of memory may be not efficient.
Because the format of the adjacency matrix is `int adjacency[nodes][degree]`, which is commanly used in regular and non-regular graphs.
This wastes memory when the number of edges that each vertex has varies greatly.

## Acknowledgment
This software is inherited from https://github.com/ryuhei-mori/graph_ASPL.git.
