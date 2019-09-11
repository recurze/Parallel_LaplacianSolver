#include "lsolver.h"

#include <cmath>
#include <random>
#include <cassert>
#include <iostream>
#include <algorithm>

double computeErrorForLaplacian(const Graph *g, const double *b, double *x) {
    int n = g->getNumVertex();
    double **L = new double*[n];
    for (int i = 0; i < n; ++i) {
        L[i] = new double[n];
    }
    g->copyLaplacianMatrix(L);

    // RMS of (Lx - b)
    double sumOfSquareError = 0;
    for (int i = 0; i < n; ++i) {
        double error = -b[i];
        for (int j = 0; j < n; ++j) {
            error += L[i][j]*x[j];
        }
        sumOfSquareError += (error * error);
        delete[] L[i];
    }
    delete[] L;
    return sqrt(sumOfSquareError/n);
}

void Lsolver::solve(const Graph *g, const double *b, double **x) {
    double *eta = NULL;
    auto beta = computeQueueOccupancyProbabilityAtStationarity(g, b, &eta);
    assert(eta != NULL);
    assert(beta > 0 and beta <= 1);

    computeCanonicalSolution(g, b, eta, beta, x);
    delete[] eta; eta = NULL;

    auto error = computeErrorForLaplacian(g, b, *x);
    std::cerr << "Beta: " << beta << "\nLaplacian Error: " << error << std::endl;
}

void Lsolver::computeJ(int n, const double *b, double *J) {
    for (int i = 0; i < n; ++i) {
        J[i] = -b[i]/b[n - 1];
    }
}

template <typename T>
T max(int n, const T *a) {
    T max_a = 0;
    for (int i = 0; i < n; ++i) {
        if (max_a < a[i]) {
            max_a = a[i];
        }
    }
    return max_a;
}

template <typename T>
T sum(int n, const T *a) {
    T sum_a = 0;
    for (int i = 0; i < n; ++i) {
        sum_a += a[i];
    }
    return sum_a;
}

template <typename T>
void addArray(int n, T *a, const T *b) {
    for (int i = 0; i < n; ++i) {
        a[i] += b[i];
    }
}

template <typename T>
void del(int n, T **P) {
    for (int i = 0; i < n; ++i) {
        delete[] P[i];
    }
    delete[] P;
}

std::mt19937 rng;
std::uniform_real_distribution<double> dist(0, 1);
inline bool trueWithProbability(double p) {
    return dist(rng) <= p;
}

void Lsolver::generateNewPackets(
        int n, int *Q, double beta, const double *J) {
    for (int i = 0; i < n - 1; ++i) {
        if (trueWithProbability(beta * J[i])) {
            ++Q[i];
        }
    }
}

void Lsolver::transmitPackets(
        int n, double **P, int *Q, int *inQ) {
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; Q[i] > 0 and j < n; ++j) {
            if (trueWithProbability(P[i][j])) {
                --Q[i];
                ++inQ[j];
            }
        }
    }
}

bool Lsolver::hasConverged(int n, int *inQ) {
    int numberOfNodesWithUnstableQueue = 0;
    for (int i = 0; i < n; ++i) {
        if (inQ[i] > 0.05*n) {
            ++numberOfNodesWithUnstableQueue;
        }
    }
    return numberOfNodesWithUnstableQueue < 0.05*n;
}

void Lsolver::updateCnt(int n, int *Q, int *cnt) {
    for (int i = 0; i < n - 1; ++i) {
        if (Q[i] > 0) {
            ++cnt[i];
        }
    }
}

void Lsolver::estimateQueueOccupancyProbability(
        int n, double **P, int *cnt, int *Q, int *inQ,
        double beta, const double *J, double T_samp, double *eta) {

    rng.seed(std::random_device{}());

    std::fill(Q, Q + n, 0);
    std::fill(cnt, cnt + n, 0);

    int T = 0;
    bool converged = false;
    bool completed = false; // completed when converged and sampled
    do {
        std::fill(inQ, inQ + n, 0);

        generateNewPackets(n, Q, beta, J);
        transmitPackets(n, P, Q, inQ);
        addArray(n, Q, inQ);

        if (not converged) {
            converged = hasConverged(n, inQ);
        } else {
            ++T;
            updateCnt(n, Q, cnt);
            completed = (T > T_samp);
        }
    } while(!completed);

    for (int i = 0; i < n; ++i) {
        eta[i] = cnt[i]/T_samp;
    }
}

template <typename T>
void initNewMemory2d(int n, int m, T*** A) {
    *A = new T*[n];
    for (int i = 0; i < n; ++i) {
        (*A)[i] = new T[m];
    }
}

double computeErrorForDCP(int n, double *eta, double **P, double beta, double *J) {
    // RMS of ((I - P)^T X eta = beta J)
    double sumOfSquareError = 0;
    for (int i = 0; i < n; ++i) {
        double error = -beta * J[i];
        for (int j = 0; j < n; ++j) {
            error += ((i == j) - P[i][j]) * eta[j];
        }
        sumOfSquareError += error * error;
    }
    return sqrt(sumOfSquareError/n);
}

double Lsolver::computeQueueOccupancyProbabilityAtStationarity(
        const Graph *g, const double *b, double **eta) {

    int n = g->getNumVertex();

    auto T_samp = 4*log(n) / (k*k*e2*e2);

    double *J = new double[n];
    computeJ(n, b, J);

    double **P = NULL;
    initNewMemory2d(n, n, &P);
    g->copyTransitionMatrix(P);

    int *cnt  = new int[n];

    // since transmission is concurrent, we need separate inbox
    int *Q = new int[n];
    int *inQ = new int[n];

    *eta = new double[n];
    double beta = 1, max_eta;
    do {
        beta /= 2;

        estimateQueueOccupancyProbability(
                n, P, cnt, Q, inQ, beta, J, T_samp, *eta);

        std::cerr << "Beta = " << beta << std::endl;
        max_eta = max(n, *eta);
    } while (max_eta > 0.75 * (1 - e1 - e2) and beta > 0);

    delete[] Q; Q = NULL;
    delete[] cnt; cnt = NULL;
    delete[] inQ; inQ = NULL;

    auto error = computeErrorForDCP(n, *eta, P, beta, J);
    std::cerr << "DCP Error: " << error << std::endl;

    del(n, P);
    delete[] J; J = NULL;

    return beta;
}

double Lsolver::computeZstar(int n, const double *eta, const double *d) {
    double zstar = 0;
    for (int i = 0; i < n; ++i) {
        zstar += eta[i]/d[i];
    }
    auto sum_d = sum(n, d);
    return -zstar*(sum_d/n);
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
        (*x)[i] = (-b[n - 1]/beta) * (eta[i]/d[i] + zstar/sum_d);
    }

    delete[] d;
    std::cerr << "Sum of x: " << sum(n, *x) << std::endl;
}
