#!/usr/bin/env python3

import sys
import numpy as np

def readInputFile(ifname):
    with open(ifname) as f:
        content = [x.strip() for x in f.readlines()]

        n, m = [int(x) for x in content[0].split()]
        A = [[0 for j in range(n)] for i in range(n)]
        for i in range(m):
            u, v, w = content[i + 1].split()
            u, v, w = int(u)-1, int(v)-1, float(w)
            if u != v:
                A[u][v] += w
                A[v][u] += w
        b = [float(x) for x in content[m + 1].split()]
        return A, b

def computeSolutionLstsq(A, b):
    L = np.array(np.diag([sum(row) for row in A]) - np.array(A))
    b = np.array(b)

    x = np.linalg.lstsq(L, b, rcond = None)[0]
    assert np.allclose(np.dot(L, x), b)
    return x

def computeSolutionJacobi(A, b):
    n = len(b)
    di = np.diag([1/sum(row) for row in A])
    A = np.array(A)
    b = np.array(b)

    x = np.zeros(n)
    while 1:
        y = di.dot(A.dot(x) + b)
        d = np.absolute(x - y)
        e = np.linalg.norm(d)
        if e < 1e-3: break
        x = y

    return x

def writeToFile(afname, x):
    n = int(afname.split('.')[0])
    with open(afname, 'w') as f:
        for i in range(n):
            print(x[i], end = " " if i < n - 1 else '', file = f)

if __name__ == "__main__":
    ifname = sys.argv[1]

    A, b = readInputFile(ifname)
    #x = computeSolutionLstsq(A, b)
    y = computeSolutionJacobi(A, b)

    afname = ifname.replace('inp', 'ans')
    writeToFile(afname, y)
