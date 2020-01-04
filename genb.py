#!/usr/bin/env python3

import random
import argparse

def generate_b(n, nsources):
    b = [0]*n
    for i in random.sample(range(n - 1), nsources):
        b[i] = 1
    b[-1] = -sum(b)
    return b

def writeToFile(b, fname):
    with open(fname, 'w') as f:
        for i in b:
            print(i, end=' ',file=f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description=
            'Generate RHS of Lx=b')

    parser.add_argument(
            'size',
            type=int,
            help='Number of vertices')

    parser.add_argument(
            '--fsources',
            type=float,
            help='Fraction of sources',
            default=0.1)


    parser.add_argument(
            '--output',
            type=str,
            help='Output file name',
            default='b$size.inp')

    args = parser.parse_args()

    n = args.size
    b = generate_b(n, int(n*args.fsources))

    ofname = args.output
    if ofname == 'b$size.inp':
        ofname = 'b' + str(n) + '.inp'

    writeToFile(b, ofname)
