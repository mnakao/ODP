#!/usr/bin/python2.7
# coding: utf-8
#=======================================================================
#
#	Create a random graph
#	by Ikki Fujiwara, National Institute of Informatics
#	2015-06-23
#
#=======================================================================
# "create-random.py" is licensed under a Creative Commons Attribution 4.0 International License.
# http://creativecommons.org/licenses/by/4.0/

author = "(random)"
email = "graphgolf@nii.ac.jp"
text1 = "A random graph provided as a baseline."

import networkx as nx
import argparse
argumentparser = argparse.ArgumentParser()
argumentparser.add_argument('nnodes', type=int)
argumentparser.add_argument('degree', type=int)
argumentparser.add_argument('-l', action='store_true')
argumentparser.add_argument('-W', required=False, type=int)
argumentparser.add_argument('-H', required=False, type=int)

def main(args):
	nnodes = args.nnodes
	degree = args.degree
        width  = 0
        heigth = 0
        if args.l:
            assert args.W != None
            assert args.H != None
            width  = args.W
            heigth = args.H
            assert nnodes == width * heigth

	assert degree < nnodes
	
	low_diam, low_aspl = lower_bound_of_diam_aspl(nnodes, degree)
	g = nx.random_regular_graph(degree, nnodes, 0)
	if nx.is_connected(g):
		hops = nx.shortest_path_length(g, weight=None)
		diam, aspl = max_avg_for_matrix(dict(hops))
	else:
		diam, aspl = float("inf"), float("inf")

        outfname = ""
        if args.l:
            print("Data Type: Grid (W x H = {} x {})".format(width, heigth))
            outfname = "w{}h{}d{}.random".format(width, heigth, degree) + ".edges"
        else:
            print("Data Type: General")
            outfname = "n{}d{}.random".format(nnodes, degree) + ".edges"

        print("Nodes        : {}".format(nnodes))
        print("Degree       : {}".format(degree))
        print("Diameter     : {}".format(diam))
        print("ASPL         : {}".format(aspl))
        print("Diameter Gap : {}".format(diam - low_diam))
        print("ASPL Gap     : {}".format(aspl - low_aspl))
        print("Output file  : {}".format(outfname))
	
        fp = open(outfname, mode='w')
        for a, b in g.edges():
            if args.l:
                line = str(a/heigth)+","+str(a%heigth)+" "+str(b/heigth)+","+str(b%heigth)+"\n"
            else:
                line = str(a)+" "+str(b)+"\n"
            fp.write(line)
            
        fp.close()
	return

def lower_bound_of_diam_aspl(nnodes, degree):
	diam = -1
	aspl = 0.0
	n = 1
	r = 1
	while True:
		tmp = n + degree * pow(degree - 1, r - 1)
		if tmp >= nnodes:
			break
		n = tmp
		aspl += r * degree * pow(degree - 1, r - 1)
		diam = r
		r += 1
	diam += 1
	aspl += diam * (nnodes - n)
	aspl /= (nnodes - 1)
	return diam, aspl

def max_avg_for_matrix(data):
	cnt = 0
	sum = 0.0
	max = 0.0

	for i in data:
		for j in data[i]:
			if i != j:
				cnt += 1
				sum += data[i][j]
				if max < data[i][j]:
					max = data[i][j]
	return max, sum / cnt

if __name__ == '__main__':
	main(argumentparser.parse_args())
