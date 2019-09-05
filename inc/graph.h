#ifndef GRAPH_H
#define GRAPH_H

class Graph {
public:
    Graph(int _n);
    Graph(int _n, double **_A);
    int getNumVertex() { return n; }
    void copyDegreeMatrix(double **_D);
    void copyLaplacianMatrix(double **_L);
    void copyAdjacencyMatrix(double **_A);
    void copyTransitionMatrix(double **_P);
private:
    int n;
    double *D;
    double **A;
    double **P;
    double **L;
};

#endif
