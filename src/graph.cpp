#include "graph.h"

Graph::Graph(int _n): n(_n) {
    D = new double[n];
    A = new double*[n];
    P = new double*[n];
    L = new double*[n];
    for (int i = 0; i < n; ++i) {
        A[i] = new double[n];
        P[i] = new double[n];
        L[i] = new double[n];
    }
}

Graph::Graph(int _n, double **_A): n(_n) {
    D = new double[n];
    A = new double*[n];
    P = new double*[n];
    L = new double*[n];
    for (int i = 0; i < n; ++i) {
        A[i] = new double[n];
        P[i] = new double[n];
        L[i] = new double[n];

        D[i] = 0;
        for (int j = 0; j < n; ++j) {
            auto x = _A[i][j];
            D[i] += x;
            A[i][j] = x;
            L[i][j] = -x;
        }
        L[i][i] = D[i];

        for (int j = 0; j < n; ++j) {
            P[i][j] = A[i][j]/D[i];
        }
    }
}

Graph::~Graph() {
    delete[] D;
    for (int i = 0; i < n; ++i) {
        delete[] A[i];
        delete[] P[i];
        delete[] L[i];
    }
    delete[] A;
    delete[] P;
    delete[] L;
}

void Graph::copyDegreeMatrix(double *_D) const {
    for (int i = 0; i < n; ++i) {
        _D[i] = D[i];
    }
}

template <typename T>
void copy2dMatrix(T **original, T **copy, int n, int m) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            copy[i][j] = original[i][j];
        }
    }
}

void Graph::copyLaplacianMatrix(double **_L) const {
    copy2dMatrix(L, _L, n, n);
}

void Graph::copyAdjacencyMatrix(double **_A) const {
    copy2dMatrix(A, _A, n, n);
}

void Graph::copyTransitionMatrix(double **_P) const {
    copy2dMatrix(P, _P, n, n);
}
