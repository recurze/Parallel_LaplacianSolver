#include "graph.h"
#include "lsolver.h"

#include <omp.h>

#include <cmath>
#include <vector>
#include <cassert>
#include <fstream>
#include <numeric>
#include <iostream>
#include <functional>

void in(const char *fname, Graph **g, std::vector<double>& b, char type);
void out(const char *fname, const std::vector<double>& x);
void computeRelError(
        Graph* g, const std::vector<double>& b, const std::vector<double>& x);

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "Usage:\n ./main <input_filename> <output_filename> 'm'/'e'\n";
        exit(0);
    }

    Graph *g = NULL;
    std::vector<double> b;
    char *ifname = argv[1];

    char type = argv[3][0];
    in(ifname, &g, b, type);

    std::vector<double> x;
    Lsolver(g, b).solve(x);

    computeRelError(g, b, x);
    delete g; g = NULL;

    char *ofname = argv[2];
    out(ofname, x);


    return 0;
}

const double EPS = 1e-6;
bool isSymmetric(const std::vector< std::vector<double> >& A) {
    int n = (int) A.size();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (fabs(A[i][j] - A[j][i]) > EPS) {
                return false;
            }
        }
    }
    return true;
}

bool isPositiveWeighted(const std::vector< std::vector<double> >& A) {
    int n = (int) A.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (A[i][j] < -EPS) {
                return false;
            }
        }
    }
    return true;
}

bool noSelfLoops(const std::vector< std::vector<double> >& A) {
    int n = (int) A.size();
    for (int i = 0; i < n; ++i) {
        if (A[i][i] > EPS) {
            return false;
        }
    }
    return true;
}

bool isConnected(const std::vector< std::vector<double> >& A) {
    int n = (int) A.size();
    std::vector<bool> visited(n, false);

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

    return true;
}

inline void checkValidGraph(const std::vector< std::vector<double> >& A) {
    assert(isSymmetric(A));
    assert(isPositiveWeighted(A));
    assert(isConnected(A));
    assert(noSelfLoops(A));
}

template <typename T>
inline T sum(const std::vector<T>& a) {
    return std::accumulate(a.begin(), a.end(), (T) 0);
}

void checkValidb(int n, const std::vector<double>& b) {
    assert(fabs(b.back()) > EPS);
    assert(fabs(sum(b)) < EPS);
}

void in(const char *fname, Graph **g, std::vector<double>& b, char type) {
    std::ifstream infile(fname);

    int n;
    infile >> n;
    assert(n > 0);

    std::vector< std::vector<double> > A(n);
    for (int i = 0; i < n; ++i) A[i].resize(n);

    if (type == 'm') {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                infile >> A[i][j];
            }
        }
    } else {
        int m;
        infile >> m;
        for (int i = 0; i < m; ++i) {
            int u, v;
            double w;
            infile >> u >> v >> w;
            if (u != v) {
                A[u - 1][v - 1] += w;
                A[v - 1][u - 1] += w;
            }
        }
    }

    checkValidGraph(A);
    *g = new Graph(A);

    b.resize(n);
    for (int i = 0; i < n; ++i) {
        infile >> b[i];
    }
    checkValidb(n, b);
}

void out(const char *fname, const std::vector<double>& x) {
    std::ofstream outfile(fname);
    for (const auto& i: x) {
        outfile << i << ' ';
    }
    outfile << '\n';
}

void computeRelError(
        Graph* g, const std::vector<double>& b, const std::vector<double>& x) {

    auto n = g->getNumVertex();
    auto D = g->getDegreeMatrix();
    auto P = g->getTransitionMatrix();

    std::vector<double> b_hat(n);

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double s = -x[i];
        for (int j = 0; j < n; ++j) {
            s += P[i][j]*x[j];
        }
        b_hat[i] = -s*D[i];
    }

    double num = 0, den = 0;
    for (int i = 0; i < n; ++i) {
        num += (b_hat[i] - b[i])*(b_hat[i] - b[i]);
        den += b[i] * b[i];
    }
    std::cerr << "Rel Error: " << num/den << '\n';
}

