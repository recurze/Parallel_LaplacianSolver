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
            P[i][j] = -L[i][j]/D[i];
        }
        P[i][i] = 0;
    }
}

void Graph::copyDegreeMatrix(double **_D) const {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            _D[i][j] = 0;
        }
        _D[i][i] = D[i];
    }
}

void Graph::copyLaplacianMatrix(double **_L) const {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            _L[i][j] = L[i][j];
        }
    }
}

void Graph::copyAdjacencyMatrix(double **_A) const {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            _A[i][j] = A[i][j];
        }
    }
}

void Graph::copyTransitionMatrix(double **_P) const {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            _P[i][j] = P[i][j];
        }
    }
}
