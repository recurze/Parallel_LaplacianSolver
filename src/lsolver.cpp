#include "lsolver.h"

#include <cmath>
#include <random>
#include <cassert>
#include <iostream>
#include <algorithm>

double computeError(const Graph *g, const double *b, double *x) {
    int n = g->getNumVertex();
    double **L = new double*[n];
    for (int i = 0; i < n; ++i) {
        L[i] = new double[n];
    }
    g->copyLaplacianMatrix(L);

    double mse = 0;
    for (int i = 0; i < n; ++i) {
        double se = -b[i];
        for (int j = 0; j < n; ++j) {
            se += L[i][j]*x[j];
        }
        mse += (se * se);
    }
    return sqrt(mse/n);
}

void Lsolver::solve(const Graph *g, const double *b, double **x) {
    double *eta = NULL;
    auto beta = computeStationaryState(g, b, &eta);
    assert(eta != NULL);
    assert(beta > 0 and beta <= 1);

    computeCanonicalSolution(g, b, eta, beta, x);
    delete[] eta; eta = NULL;

    auto error = computeError(g, b, *x);
    std::cerr << "Beta: " << beta << "\nError: " << error << std::endl;
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
void del(int n, T **P) {
    for (int i = 0; i < n; ++i) {
        delete[] P[i];
    }
    delete[] P;
}

template <typename T>
void copy1d(int n, T *orig, T *copy) {
    for (int i = 0; i < n; ++i) {
        copy[i] = orig[i];
    }
}

std::mt19937 rng;
std::uniform_real_distribution<double> dist(0, 1);

inline bool generatePacket(double p) {
    return dist(rng) <= p;
}

// If rand lands in [prefixPi[i], prefixPi[i + 1]) which is probability Pi[i],
int pickNeighbor(int n, double *prefixPi) {
    auto x = std::upper_bound(prefixPi, prefixPi + n, dist(rng)) - prefixPi;
    return x;
}

void Lsolver::generateNewPackets(
        int n, int *Q, double beta, const double *J) {
    for (int i = 0; i < n - 1; ++i) {
        Q[i] += generatePacket(beta * J[i]);
    }
}

void Lsolver::transmitToRandomNeighbor(
        int n, double *prefixPi, int *Q, int qid) {
    --Q[qid];
    ++Q[pickNeighbor(n, prefixPi)];
}

void Lsolver::transmitPackets(
        int n, double **prefixP, int *oldQ, int *newQ) {
    for (int i = 0; i < n - 1; ++i) {
        if (oldQ[i] > 0) {
            transmitToRandomNeighbor(n, prefixP[i], newQ, i);
        }
    }
}

bool Lsolver::hasConverged(int n, int *oldQ, int *newQ) {
    int changed = 0;
    for (int i = 0; i < n; ++i) {
        if (oldQ[i] != newQ[i]) {
            ++changed;
        }
    }
    return changed < 0.1*n;
}

void Lsolver::updateCnt(int n, int *Q, int *cnt) {
    for (int i = 0; i < n - 1; ++i) {
        if (Q[i] > 0) {
            ++cnt[i];
        }
    }
}

void Lsolver::estimateQueueOccupancyProbability(
        int n, double **prefixP, int *cnt, int *oldQ, int *newQ,
        double beta, const double *J, double T_samp, double *eta) {

    rng.seed(std::random_device{}());

    std::fill(cnt, cnt + n, 0);
    std::fill(oldQ, oldQ + n, 0);
    std::fill(newQ, newQ + n, 0);

    int T = 0;
    bool converged = false;
    bool completed = false;
    do {
        generateNewPackets(n, oldQ, beta, J);
        copy1d(n, oldQ, newQ);

        transmitPackets(n, prefixP, oldQ, newQ);
        copy1d(n, newQ, oldQ);

        if (not converged) {
            converged = hasConverged(n, oldQ, newQ);
        } else {
            ++T;
            updateCnt(n, oldQ, cnt);
            completed = (T > T_samp);
        }
    } while(!completed);

    for (int i = 0; i < n; ++i) {
        eta[i] = cnt[i]/T_samp;
    }
}

template <typename T>
void init2d(int n, int m, T*** A) {
    *A = new T*[n];
    for (int i = 0; i < n; ++i) {
        (*A)[i] = new T[m];
    }
}

template <typename T>
void makePrefixSum(int n, T* A) {
    for (int i = 1; i < n; ++i) {
        A[i] += A[i - 1];
    }
}

void Lsolver::computePrefixP(int n, const Graph *g, double **prefixP) {
    g->copyTransitionMatrix(prefixP);
    // check pickNeighbor to see it work in log(n) thanks to prefixSum of P
    for (int i = 0; i < n; ++i) {
        makePrefixSum(n, prefixP[i]);
    }
}

double checkStationarity(int n, double **P, double *pi) {
    double mse = 0;
    for (int i = 0; i < n; ++i) {
        double se = -pi[i];
        for (int j = 0; j < n; ++j) {
            se += (P[j][i] * pi[j]);
        }
        mse += se * se;
    }
    return sqrt(mse/n);
}

double Lsolver::computeStationaryState(
        const Graph *g, const double *b, double **eta) {

    int n = g->getNumVertex();

    auto T_samp = 4*log(n) / (k*k*e2*e2);

    double *J = new double[n];
    computeJ(n, b, J);

    double **prefixP = NULL;
    init2d(n, n, &prefixP);
    computePrefixP(n, g, prefixP);

    int *cnt = new int[n];
    int *oldQ = new int[n];
    int *newQ = new int[n];

    *eta = new double[n];
    double beta = 1, max_eta;
    do {
        beta /= 2;

        estimateQueueOccupancyProbability(
                n, prefixP, cnt, oldQ, newQ, beta, J, T_samp, *eta);

        max_eta = max(n, *eta);
    } while (max_eta > 0.75 * (1 - e1 - e2) and beta > 0);

    delete[] J; J = NULL;
    delete[] cnt; cnt = NULL;
    delete[] oldQ; oldQ = NULL;
    delete[] newQ; newQ = NULL;

    g->copyTransitionMatrix(prefixP);
    auto error = checkStationarity(n, prefixP, *eta);
    del(n, prefixP);

    std::cerr << "Stationarity: " << error << std::endl;
    return beta;
}

double Lsolver::computeZstar(int n, const double *eta, const double *d) {
    double zstar = 0;
    for (int i = 0; i < n; ++i) {
        zstar += eta[i]/d[i];
    }
    return -zstar;
}

void Lsolver::computeCanonicalSolution(
        const Graph *g, const double *b,
        double *eta, double beta, double **x) {

    int n = g->getNumVertex();

    double* d = new double[n];
    g->copyDegreeMatrix(d);

    auto zstar = computeZstar(n, eta, d);

    auto sum_d = sum(n, d);

    *x = new double[n];
    for (int i = 0; i < n; ++i) {
        (*x)[i] = (-b[n - 1]/beta) * (eta[i]/d[i] + zstar*(d[i]/sum_d));
    }
    std::cerr << "Sum of x: " << sum(n, *x) << std::endl;
}
