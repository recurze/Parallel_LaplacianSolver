#include "lsolver.h"

#include <cmath>

void Lsolver::solve(const Graph *g, const double *b, double *x) {
    double *eta = NULL;
    auto beta = computeStationaryState(g, b, eta);

    computeCanonicalSolution(g, b, eta, beta, x);

    delete[] eta; eta = NULL;
}

void Lsolver::computeJ(int n, const double *b, double *J) {
    for (int i = 0; i < n; ++i) {
        J[i] = -b[i]/b[n - 1];
    }
}

double max(int n, double *eta) {
    double max_eta = 0;
    for (int i = 0; i < n; ++i) {
        if (max_eta < eta[i]) {
            max_eta = eta[i];
        }
    }
    return max_eta;
}

template <typename T>
void init2dMatrix(T **P, int n, int m) {
    P = new T*[n];
    for (int i = 0; i < n; ++i) {
        P[i] = new T[m];
    }
}

template <typename T>
void del(T **P, int n) {
    for (int i = 0; i < n; ++i) {
        delete[] P[i];
    }
    delete[] P;
}

// TODO: This is the main parallel function!
void Lsolver::estimateQueueOccupancyProbability(
        double **P, double beta,
        const double *J, double T_samp, double *eta) {
}

double Lsolver::computeStationaryState(
        const Graph *g, const double *b, double *eta) {
    int n = g->getNumVertex();

    auto T_samp = 4*log(n) / (k*k*e2*e2);

    double *J = new double[n];
    computeJ(n, b, J);

    double **P = NULL;
    init2dMatrix(P, n, n);
    g->copyTransitionMatrix(P);

    double beta = 1;
    eta = new double[n];
    do {
        beta /= 2;
        estimateQueueOccupancyProbability(P, beta, J, T_samp, eta);
    } while (max(n, eta) < 0.75 * (1 - e1 - e2));

    del(P, n);
    delete[] J; J = NULL;

    return beta;
}

double Lsolver::computeZstar(int n, const double *eta, const double *d) {
    double zstar = 0;
    for (int i = 0; i < n; ++i) {
        zstar += eta[i]/d[i];
    }
    return -zstar;
}


// TODO: This can also be parallely after computing z*
void Lsolver::computeCanonicalSolution(
        const Graph *g, const double *b,
        double *eta, double beta, double *x) {
    int n = g->getNumVertex();

    double d = new double[n];
    g->copyDegreeMatrix(d);

    auto zstar = computeZstar(n, eta, d);

    x = new double[n];
}
