#ifndef GRAPH_H
#define GRAPH_H

#include <vector>

class Graph {
public:
    Graph(const std::vector< std::vector<double> >& _A);

    int getNumVertex() const { return n; }
    std::vector<double> getDegreeMatrix() const { return D; }
    std::vector< std::vector<double> > getTransitionMatrix() const { return P; }
private:
    int n;
    std::vector<double> D;
    std::vector< std::vector<double> > P;
};

#endif
