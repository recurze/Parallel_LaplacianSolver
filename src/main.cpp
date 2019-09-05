#include "graph.h"
#include "lsolver.h"

#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>

void in(const char *fname, Graph *g, double *b);
void out(const char *fname, int n, double *x, double e, double t);
double computeError(int n, double **L, const double *x, const double *b);

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage:\n ./main <input_filename> <output_filename>\n";
        exit(0);
    }

    Graph *g = NULL;
    double *b = NULL;
    char *ifname = argv[0];
    in(ifname, g, b);

    auto ts = std::chrono::high_resolution_clock::now();
    double *x = NULL;
    Lsolver().solve(g, b, x);
    auto te = std::chrono::high_resolution_clock::now();


    auto n = g->getNumVertex();

    double **L = new double*[n];
    for (int i = 0; i < n; ++i) {
        L[i] = new double[n];
    }
    g->copyLaplacianMatrix(L);
    auto error = computeError(n, L, x, b);

    auto time = std::chrono::duration_cast<std::chrono::microseconds>(te - ts).count();

    char *ofname = argv[1];
    out(ofname, n, x, error, time);

    return 0;
}

void in(const char *fname, Graph *g, double *b) {
    std::ifstream infile(fname);

    int numVertex;
    infile >> numVertex;

    b = new double[numVertex];

    double **A = new double*[numVertex];
    for (int i = 0; i < numVertex; ++i) {
        A[i] = new double[numVertex];
        for (int j = 0; j < numVertex; ++j) {
            infile >> A[i][j];
        }
        infile >> b[i];
    }

    g = new Graph(numVertex, A);
}

double computeError(int n, const double **L, const double *x, const double *b) {
    double mse = 0;
    for (int i = 0; i < n; ++i) {
        double e = -b[i];
        for (int j = 0; j < n; ++j) {
            e += L[i][j] * x[j];
        }
        mse += e*e;
    }
    mse /= n;
    return sqrt(mse);
}

void out(const char *fname, int n, double *x, double e, double t) {
    std::ofstream outfile(fname);
    outfile << "Error: " << e << '\n';
    outfile << "Time taken: " << t << '\n';

    for (int i = 0; i < n; ++i) {
        outfile << x[i] << (i == n - 1 ? '\n' : ' ');
    }
}

