#!/usr/bin/env python3

'''
Purpose:
    Generate undirected positive weighted connected graph of given size

Command line argument:
    n: integer denoting number of nodes in the graph

Console Output:
    If not run properly, instructions on usage; else none

File output:
    filename: n.inp; replace n with integer
    Contents: n + 2 lines; format as specified in the README
'''

import sys
import random
import numpy as np
from math import sqrt

minw, maxw = 0.01, 100
def generate_undirectedWeightedConnectedGraph(n, m):
    assert m >= n - 1 and m <= n*(n - 1)/2

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

def writeToFile(n, A, b):
    ifname = str(n) + '.inp'
    with open(ifname, 'w') as f:
        print(n, file=f)
        for i in range(n):
            for j in range(n):
                print(A[i][j], end = " " if j < n - 1 else '', file = f)
            print("", file = f)
        for i in range(n):
            print(b[i], end = " " if i < n - 1 else '', file = f)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Example Usage: ./testgen.py 10")
        exit(0)
    n = int(sys.argv[1])
    m = random.randint(n - 1, n*(n - 1)//2)
    A = generate_undirectedWeightedConnectedGraph(n, m)
    b = generate_randomb(n)

    writeToFile(n, A, b)

