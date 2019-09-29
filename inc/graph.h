#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <iostream>

class Graph {
public:
    Graph() {}
    Graph(int _n): n(_n) {
        init();
    }

    void setNumVertex(int _n) {
        n = _n;
        init();
    }

    int getNumVertex() const {
        return n;
    }

    void addEdge(int u, int v, double w) {
        deg[u] += w;
        adj[u].push_back({v, w});

        deg[v] += w;
        adj[v].push_back({u, w});
    }

    const std::vector< std::pair<int, double> >& getNeighbors(int u) const {
        return adj[u];
    }

    const std::vector<double>& getDegreeMatrix() const {
        return deg;
    }
private:
    int n;

    std::vector<double> deg;
    std::vector< std::vector< std::pair<int, double> > > adj;

    void init() {
        deg.resize(n, 0);
        adj.resize(n);
    }
};

#endif
