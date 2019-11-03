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

minw, maxw = 1, 100
def generate_undirectedWeightedConnectedGraph(n, m):
    assert m >= n - 1 and m <= n*(n - 1)/2

    nodes = list(range(0, n))
    A = [{} for j in range(n)]

    def add_random_edges(m):
        n_edges = 0
        for _ in range(m):
            i, j = random.sample(nodes, 2)
            if j not in A[i] and i not in A[j]:
                w = random.uniform(minw, maxw)
                A[i][j] = w
                n_edges += 1
        return n_edges

    def mst():
        S, T = set(nodes), set()

        curr = random.sample(S, 1).pop()

        T.add(curr)
        S.remove(curr)
        while S:
            new = random.sample(nodes, 1).pop()
            if new not in T:
                w = random.uniform(minw, maxw)
                A[curr][new] = w
                T.add(new)
                S.remove(new)
            curr = new

    mst()
    n_edges = add_random_edges(m - n + 1)
    return A, n_edges + n - 1

minb, maxb = 0, 1
def generate_randomb(n):
    b = [random.randint(minb, maxb) for _ in range(n - 1)]
    b.append(-sum(b))
    return b

def writeToFile(n, m, A, b):
    ifname = str(n) + '.inp'
    with open(ifname, 'w') as f:
        print(n, m, file=f)
        for i in range(n):
            for j in sorted(A[i]):
                print(i + 1, j + 1, A[i][j], file = f)
        for i in range(n):
            print(b[i], end = " ", file = f)

if __name__ == "__main__":
    n = int(sys.argv[1])
    if len(sys.argv) == 3:
        m = int(n * float(sys.argv[2]))
    else:
        m = random.randint(n - 1, n*(n - 1)//2)
    A, m = generate_undirectedWeightedConnectedGraph(n, m)
    b = generate_randomb(n)

    writeToFile(n, m, A, b)

