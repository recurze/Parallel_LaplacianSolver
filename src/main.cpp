#include "graph.h"
#include "lsolver.h"

#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>
#include <functional>

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
    Lsolver(0.1, 0.1, 0.1).solve(g, b, &x);

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

bool isConnected(int n, double **A) {
    bool *visited = new bool[n];
    std::fill(visited, visited + n, false);

    std::function<void(int)> dfs = [&](int u) {
        visited[u] = true;
        for (int v = 0; v < n; ++v) {
            if (A[u][v] > EPS and not visited[v]) {
                dfs(v);
            }
        }
    };

    dfs(0);
    for (int i = 0; i < n; ++i) {
        if (not visited[i]) {
            return false;
        }
    }

    delete[] visited;
    return true;
}

inline void checkValidGraph(int n, double **A) {
    assert(isSymmetric(n, A));
    assert(isPositiveWeighted(n, A));
    assert(isConnected(n, A));
    assert(noSelfLoops(n, A));
}

template <typename T>
T sum(int n, const T *a) {
    T sum_a = 0;
    for (int i = 0; i < n; ++i) {
        sum_a += a[i];
    }
    return sum_a;
}

void checkValidb(int n, const double *b) {
    assert(fabs(b[n - 1]) > EPS);
    assert(fabs(sum(n, b)) < EPS);
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

    for (int i = 0; i < n; ++i) {
        delete[] A[i];
    }
    delete[] A;

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
