#!/usr/bin/env python3

'''
Generates 1 files
    n.inp
    where n is the size of the graph, fed as command line argument

    n.inp contains the graph and b according to the format specified in README
'''

import sys
import random
import numpy as np

minw, maxw = 0.1, 100.0
def generate_undirectedWeightedConnectedGraph(n, m):
    assert m >= n - 1 and m <= n*(n - 1)/2
    msg = "Generating a random undirected weighted connected graph with {} vertices and < {} edges".format(n, m + 1)
    print(msg)

    nodes = list(range(0, n))
    A = [[0 for i in range(n)] for j in range(n)]

    def add_random_edges(m):
        for _ in range(m):
            i, j = random.sample(nodes, 2)
            if A[i][j] == 0:
                w = random.uniform(minw, maxw)
                A[i][j] = A[j][i] = w

    def mst():
        S, T = set(nodes), set()

        curr = random.sample(S, 1).pop()

        T.add(curr)
        S.remove(curr)
        while S:
            new = random.sample(nodes, 1).pop()
            if new not in T:
                w = random.uniform(minw, maxw)
                A[curr][new] = A[new][curr] = w
                T.add(new)
                S.remove(new)
            curr = new

    mst()
    add_random_edges(m - n + 1)
    return A

minb, maxb = 0, 100
def generate_randomb(n):
    b = [random.uniform(minb, maxb) for _ in range(n - 1)]
    b.append(-sum(b))
    return b

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage:\n ./testgen.py n")
        exit(0)

    n = int(sys.argv[1])
    m = random.randint(n - 1, n*(n - 1)/2)
    A = generate_undirectedWeightedConnectedGraph(n, m)
    b = generate_randomb(n)

    ifname = str(n) + '.inp'
    with open(ifname, 'w') as f:
        print(n, file=f)
        for i in range(n):
            for x in A[i]:
                print(x, end = " ", file = f)
            print("", file = f)
        for x in b:
            print(x, end = " ", file = f)
