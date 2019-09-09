#include "graph.h"
#include "lsolver.h"

#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>

void in(const char *fname, Graph **g, double **b);
void out(const char *fname, int n, double *x);

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage:\n ./main <input_filename> <output_filename>\n";
        exit(0);
    }

    Graph *g = NULL;
    double *b = NULL;
    char *ifname = argv[1];
    in(ifname, &g, &b);

    assert(g != NULL);
    assert(b != NULL);

    double *x = NULL;
    Lsolver().solve(g, b, &x);

    assert(x != NULL);

    delete[] b; b = NULL;

    int n = g->getNumVertex();
    delete g; g = NULL;

    char *ofname = argv[2];
    out(ofname, n, x);

    delete[] x; x = NULL;

    return 0;
}

const double EPS = 1e-6;

bool isSymmetric(int n, double **A) {
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (fabs(A[i][j] - A[j][i]) > EPS) {
                return false;
            }
        }
    }
    return true;
}

bool isPositiveWeighted(int n, double **A) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (A[i][j] < -EPS) {
                return false;
            }
        }
    }
    return true;
}

bool noSelfLoops(int n, double **A) {
    for (int i = 0; i < n; ++i) {
        if (A[i][i] > EPS) {
            return false;
        }
    }
    return true;
}

inline void checkValidGraph(int n, double **A) {
    assert(isSymmetric(n, A));
    assert(isPositiveWeighted(n, A));
    assert(noSelfLoops(n, A));
}

void checkValidb(int n, double *b) {
    double s = 0;
    for (int i = 0; i < n; ++i) {
        s += b[i];
    }
    assert(fabs(s) < EPS);
    assert(fabs(b[n - 1]) > EPS);
}

void in(const char *fname, Graph **g, double **b) {
    std::ifstream infile(fname);

    int n;
    infile >> n;
    assert(n > 0);

    double **A = new double*[n];
    for (int i = 0; i < n; ++i) {
        A[i] = new double[n];
        for (int j = 0; j < n; ++j) {
            infile >> A[i][j];
        }
    }
    checkValidGraph(n, A);
    *g = new Graph(n, A);

    *b = new double[n];
    for (int i = 0; i < n; ++i) {
        infile >> (*b)[i];
    }
    checkValidb(n, *b);
}

void out(const char *fname, int n, double *x) {
    std::ofstream outfile(fname);
    for (int i = 0; i < n; ++i) {
        outfile << x[i] << (i == n - 1 ? '\n' : ' ');
    }
}
