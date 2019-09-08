#include "graph.h"
#include "lsolver.h"

#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>

void in(const char *fname, Graph *g, double *b);
void out(const char *fname, int n, double *x);

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage:\n ./main <input_filename> <output_filename>\n";
        exit(0);
    }

    Graph *g = NULL;
    double *b = NULL;
    char *ifname = argv[1];
    in(ifname, g, b);

    double *x = NULL;
    Lsolver().solve(g, b, x);

    delete[] b; b = NULL;

    int n = g->getNumVertex();
    delete g; g = NULL;

    char *ofname = argv[2];
    out(ofname, n, x);

    delete[] x; x = NULL;

    return 0;
}

void in(const char *fname, Graph *g, double *b) {
    std::ifstream infile(fname);

    int n;
    infile >> n;

    double **A = new double*[n];
    for (int i = 0; i < n; ++i) {
        A[i] = new double[n];
        for (int j = 0; j < n; ++j) {
            infile >> A[i][j];
        }
    }
    g = new Graph(n, A);

    b = new double[n];
    for (int i = 0; i < n; ++i) {
        infile >> b[i];
    }
}

void out(const char *fname, int n, double *x) {
    std::ofstream outfile(fname);
    for (int i = 0; i < n; ++i) {
        outfile << x[i] << (i == n - 1 ? '\n' : ' ');
    }
}
