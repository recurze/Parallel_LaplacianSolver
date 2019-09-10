#ifndef GRAPH_H
#define GRAPH_H

class Graph {
public:
    Graph(int _n);
    Graph(int _n, double **_A);
    ~Graph();

    int getNumVertex() const { return n; }
    void copyDegreeMatrix(double *_D) const;
    void copyLaplacianMatrix(double **_L) const;
    void copyAdjacencyMatrix(double **_A) const;
    void copyTransitionMatrix(double **_P) const;
private:
    int n;
    double *D;
    double **A;
    double **P;
    double **L;

    void init();
    void computeL();
    void computeP();
    void computeD();
};

#endif
