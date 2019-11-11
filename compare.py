#!/usr/bin/env python3

import sys
import numpy as np

def readOutputFile(ofname):
    with open(ofname) as f:
        content = [x.strip() for x in f.readlines()]

        x = np.array([float(i) for i in content[0].split()])
        return x

if __name__ == "__main__":
    ofname = sys.argv[1]
    afname = sys.argv[2]

    x = readOutputFile(afname)
    x_hat = readOutputFile(ofname)

    x -= np.min(x)
    x_hat -= np.min(x_hat)

    diff = np.absolute(x - x_hat)
    err = np.linalg.norm(diff)/np.linalg.norm(x)
    print(err, file=sys.stderr)
