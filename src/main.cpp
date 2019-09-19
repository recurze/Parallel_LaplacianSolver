#include "graph.h"
#include "lsolver.h"

#include <cmath>
#include <vector>
#include <cassert>
#include <fstream>
#include <numeric>
#include <iostream>
#include <functional>

void in(const char *fname, Graph **g, std::vector<double>& b);
void out(const char *fname, const std::vector<double>& x);

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage:\n ./main <input_filename> <output_filename>\n";
        exit(0);
    }

    Graph *g = NULL;
    std::vector<double> b;
    char *ifname = argv[1];

    in(ifname, &g, b);

    std::vector<double> x;
    Lsolver(g, b).solve(x);

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

void in(const char *fname, Graph **g, std::vector<double>& b) {
    std::ifstream infile(fname);

    int n;
    infile >> n;
    assert(n > 0);

    std::vector< std::vector<double> > A(n);
    for (int i = 0; i < n; ++i) {
        A[i].resize(n);
        for (int j = 0; j < n; ++j) {
            infile >> A[i][j];
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
