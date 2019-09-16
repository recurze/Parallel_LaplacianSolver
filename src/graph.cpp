#include "graph.h"

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

void Graph::init() {
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

Graph::Graph(int _n): n(_n) {
    init();
}

template <typename T>
T sum(int n, const T *a) {
    T sum_a = 0;
    for (int i = 0; i < n; ++i) {
        sum_a += a[i];
    }
    return sum_a;
}

void Graph::computeD() {
    for (int i = 0; i < n; ++i) {
        D[i] = sum(n, A[i]);
    }
}

void Graph::computeP() {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            P[i][j] = A[i][j]/D[i];
        }
    }
}

void Graph::computeL() {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            L[i][j] = -A[i][j];
        }
        L[i][i] = D[i];
    }
}

Graph::Graph(int _n, double **_A): n(_n) {
    init();
    copy2d(_A, A, n, n);
    computeD();
    computeP();
    computeL();
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
