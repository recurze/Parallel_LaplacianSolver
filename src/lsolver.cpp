#include "lsolver.h"

#include <cmath>
#include <random>
#include <algorithm>

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

template <typename T>
T max(int n, T *a) {
    T max_a = 0;
    for (int i = 0; i < n; ++i) {
        if (max_a < a[i]) {
            max_a = a[i];
        }
    }
    return max_a;
}

template <typename T>
T sum(int n, T *a) {
    T sum_a = 0;
    for (int i = 0; i < n; ++i) {
        sum_a += a[i];
    }
    return sum_a;
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

template <typename T>
void copy(int n, T *orig, T *copy) {
    for (int i = 0; i < n; ++i) {
        copy[i] = orig[i];
    }
}

std::mt19937 rng;
std::uniform_real_distribution<double> dist(0, 1);

bool generatePacket(double p) {
    return dist(rng) <= p;
}

int pickNeighbor(int n, double *P) {
    double *prefix_P = new double[n];

    prefix_P[0] = P[0];
    for (int i = 1; i < n; ++i) {
        prefix_P[i] = prefix_P[i - 1] + P[i];
    }
    auto x = std::upper_bound(prefix_P, prefix_P + n, dist(rng)) - prefix_P;

    delete[] prefix_P;
    return x;
}

// TODO: This is the main parallel function!
void Lsolver::estimateQueueOccupancyProbability(
        int n, double **P, double beta,
        const double *J, double T_samp, double *eta) {


            rng.seed(std::random_device{}());
    int * Q = new int[n];
    int *nQ = new int[n];

    bool converged = false;
    bool completed = false;
    do {
        for (int i = 0; i < n; ++i) {
            Q[i] += generatePacket(beta * J[i]);
        }

        copy(n, Q, nQ);

        for (int i = 0; i < n; ++i) {
            if (Q[i] > 0) {
                auto v = pickNeighbor(n, P[i]);
                --nQ[i];
                ++nQ[v];
            }
        }

        if (not converged) {
            int c = 0;
            for (int i = 0; i < n; ++i) {
                c += (Q[i] != nQ[i]);
            }
            converged = (c < 0.1*n);
        }
        copy(n, nQ, Q);

        if (converged) {
            T_samp -= 1;
            for (int i = 0; i < n; ++i) {
                eta[i] += (Q[i] > 0);
            }
            completed = (T_samp < 0);
        }
    } while (!completed);

    delete[]  Q;
    delete[] nQ;

    for (int i = 0; i < n; ++i) {
        eta[i] /= T_samp;
    }
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
        estimateQueueOccupancyProbability(n, P, beta, J, T_samp, eta);
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


// Parallelize last for loop
void Lsolver::computeCanonicalSolution(
        const Graph *g, const double *b,
        double *eta, double beta, double *x) {
    int n = g->getNumVertex();

    double* d = new double[n];
    g->copyDegreeMatrix(d);

    auto zstar = computeZstar(n, eta, d);

    auto sum_d = sum(n, d);

    x = new double[n];
    for (int i = 0; i < n; ++i) {
        x[i] = (-b[n - 1]/beta) * (eta[i]/d[i] + zstar*(d[i]/sum_d));
    }
}
