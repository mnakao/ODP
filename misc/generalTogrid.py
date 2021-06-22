#!/usr/bin/python2.7
# coding: utf-8
import argparse
import os
argumentparser = argparse.ArgumentParser()
argumentparser.add_argument("edges_path", help="Input edgelist file path")
argumentparser.add_argument('-W', required=True, type=int)
argumentparser.add_argument('-H', required=True, type=int)
args = argumentparser.parse_args()

def main():
    infile = os.path.normpath(args.edges_path)
    assert os.path.isfile(infile)
    lines = sum(1 for line in open(infile))
    width  = args.W
    height = args.H
#    print "Lines  : ", lines
#    print "Width  : ", width
#    print "Height : ", height

    f = open(args.edges_path, "r")
    for line in f:
        data = line.split()
        print("{},{} {},{}".format(int(data[0])/height, int(data[0])%height, int(data[1])/height, int(data[1])%height))
              
    f.close()

if __name__ == '__main__':
	main()
