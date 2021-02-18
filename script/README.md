## Create random general and grid graphs
```
create_random.py [-l] [-W W] [-H H] nodes degree
```

* -l : Create a grid graph
* -W : Width for a Grid graph
* -H : Height for a Grid graph
* nodes : Number of nodes
* degree : Degree in a graph

### Example for creating a general graph with 10 nodes and 3 degree
```
$ python ./create_random.py 10 3
Data Type: General
Nodes        : 10
Degree       : 3
Diameter     : 3
ASPL         : 1.88888888889
Diameter Gap : 1
ASPL Gap     : 0.222222222222
Output file  : n10d3.random.edges
```

### Example for	creating a general graph with 4x5 nodes and 3 degree
```
$ python ./create_random.py -l -W 4 -H 5 20 3
Data Type: Grid (W x H = 4 x 5)
Nodes        : 20
Degree       : 3
Diameter     : 5
ASPL         : 2.70526315789
Diameter Gap : 2
ASPL Gap     : 0.336842105263
Output file  : w4h5d3.random.edges
```

## Draw general and grid graphs
```
draw_graph.py [-n] [-s S] edge
```

* -n : Add number label
* -s : Image size
* edge : Edge file

## Calculate ASPL and diameter for a general graph
```
verfy_general.py edge
```

* edge : Edge file

### Example
```
$ python ./verfy_general.py n10d3.random.edges
Nodes = 10, Degrees = 3
Diameter     = 3
Diameter Gap = 1 (3 - 2)
ASPL         = 1.8888888889 (85/45)
ASPL Gap     = 0.2222222222 (1.8888888889 - 1.6666666667)
```

