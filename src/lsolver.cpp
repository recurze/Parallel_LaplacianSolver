#include "lsolver.h"

#include <cmath>

void Lsolver::solve(const Graph *g, const double *b, double *x) {
    auto t_hit = computeT_hit(g);

    double *eta = NULL;
    auto beta = computeStationaryState(g, b, t_hit, eta);

    computeCanonicalSolution(g, b, eta, beta, x);

    delete[] eta; eta = NULL;
}

// TODO: Figure out how to compute t_hit
double Lsolver::computeT_hit(const Graph *g) {
    return (double) g->getNumVertex();
}

void computeJ(int n, const double *b, double *J) {
    double b_sink = 0;
    for (int i = 0; i < n; ++i) {
        J[i] = -b[i];
        if (b[i] < 0) b_sink = b[i];
    }
    for (int i = 0; i < n; ++i) {
        J[i] /= b_sink;
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
void del(T **P, int n, int m) {
    for (int i = 0; i < n; ++i) {
        delete[] P[i];
    }
    delete[] P;
}

// TODO: This is the main parallel function!
void Lsolver::estimateQueueOccupancyProbability(
        double **P, const double beta, const double *J,
        const double T_mix, const double T_samp, double *eta) {
}

double Lsolver::computeStationaryState(
        const Graph *g, const double *b, double t_hit, double *eta) {
    int n = g->getNumVertex();

    auto T_mix = 64 * t_hit * log(1.0/e1);
    auto T_samp = 4 * log(n)/(e2*e2); // TODO: Do we need kappa?

    double *J = new double[n];
    computeJ(n, b, J);

    double **P = NULL;
    init2dMatrix(P, n, n);
    g->copyTransitionMatrix(P);

    double beta = 1;
    eta = new double[n];
    do {
        beta /= 2;
        estimateQueueOccupancyProbability(P, beta, J, T_mix, T_samp, eta);
    } while (max(n, eta) < 0.75 * (1 - e1 - e2));

    del(P, n, n);
    delete[] J; J = NULL;

    return beta;
}

// TODO
void Lsolver::computeCanonicalSolution(
        const Graph *g, const double *b, double *eta, double beta, double *x) {
    x = new double[g->getNumVertex()];
}
