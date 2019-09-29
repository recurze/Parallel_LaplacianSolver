#include "graph.h"
#include "lsolver.h"

#include <cmath>
#include <vector>
#include <cassert>
#include <fstream>
#include <numeric>
#include <iostream>
#include <functional>

void in(const char *fname, Graph& g, std::vector<double>& b);
void out(const char *fname, const std::vector<double>& x);

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage:\n ./main <input_filename> <output_filename>\n";
        exit(0);
    }

    Graph g;
    std::vector<double> b;
    char *ifname = argv[1];

    in(ifname, g, b);

    auto x = Lsolver(g, b).solve();

    char *ofname = argv[2];
    out(ofname, x);

    return 0;
}

const double EPS = 1e-6;

bool isConnected(const Graph& g) {
    auto n = g.getNumVertex();
    std::vector<bool> visited(n, false);

    std::function<void(int)> dfs = [&](int u) {
        visited[u] = true;
        for (auto& neigh: g.getNeighbors(u)) {
            auto v = neigh.first;
            if (not visited[v]) {
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

void checkValidb(const std::vector<double>& b) {
    auto b_sink = b.back();
    assert(fabs(b_sink) > EPS);

    auto sum_b = std::accumulate(b.begin(), b.end(), 0.0);
    assert(fabs(sum_b) < EPS);
}

void in(const char *ifname, Graph& g, std::vector<double>& b) {
    std::ifstream infile(ifname);

    int n;
    int m;
    infile >> n >> m;

    g.setNumVertex(n);
    for (int i = 0; i < m; ++i) {
        int u;
        int v;
        double w;

        infile >> u >> v >> w;

        assert(w > EPS);
        assert(u != v);

        g.addEdge(u - 1, v - 1, w);
    }
    assert(isConnected(g));

    b.resize(n);
    for (int i = 0; i < n; ++i) {
        infile >> b[i];
    }
    checkValidb(b);
}

void out(const char *fname, const std::vector<double>& x) {
    std::ofstream outfile(fname);
    for (const auto& i: x) {
        outfile << i << ' ';
    }
    outfile << '\n';
}
