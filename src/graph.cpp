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
    for (int i = 0; i < n; ++i) {
        delete[] A[i];
        delete[] P[i];
        delete[] L[i];
    }
    delete[] D;
    delete[] A;
    delete[] P;
    delete[] L;
}

template <typename T>
void copy1d(T *orig, T *copy, int n) {
    for (int i = 0; i < n; ++i) {
        copy[i] = orig[i];
    }
}

template <typename T>
void copy2d(T **orig, T **copy, int n, int m) {
    for (int i = 0; i < n; ++i) {
        copy1d(orig[i], copy[i], m);
    }
}

void Graph::copyDegreeMatrix(double *_D) const {
    copy1d(D, _D, n);
}

void Graph::copyLaplacianMatrix(double **_L) const {
    copy2d(L, _L, n, n);
}

void Graph::copyAdjacencyMatrix(double **_A) const {
    copy2d(A, _A, n, n);
}

void Graph::copyTransitionMatrix(double **_P) const {
    copy2d(P, _P, n, n);
}
