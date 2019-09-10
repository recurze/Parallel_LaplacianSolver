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

    // RMS of (Lx - b)
    double sumOfSquareError = 0;
    for (int i = 0; i < n; ++i) {
        double error = -b[i];
        for (int j = 0; j < n; ++j) {
            error += L[i][j]*x[j];
        }
        sumOfSquareError += (error * error);
    }
    return sqrt(sumOfSquareError/n);
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
void copy1d(T *orig, T *copy, int n) {
    for (int i = 0; i < n; ++i) {
        copy[i] = orig[i];
    }
}


std::mt19937 rng;
std::uniform_real_distribution<double> dist(0, 1);
inline bool generatePacketWithProbability(double p) {
    return dist(rng) <= p;
}

// Random lands in [prefixPi[i], prefixPi[i + 1]) with probability Pi[i]
int pickRandomNeighbor(int n, double *prefixPi) {
    auto x = std::upper_bound(prefixPi, prefixPi + n, dist(rng)) - prefixPi;
    return x;
}

void Lsolver::generateNewPackets(
        int n, int *Q, double beta, const double *J) {
    for (int i = 0; i < n - 1; ++i) {
        Q[i] += generatePacketWithProbability(beta * J[i]);
    }
}

void Lsolver::transmitToRandomNeighbor(
        int n, double *prefixPi, int *Q, int qid) {
    --Q[qid];
    ++Q[pickRandomNeighbor(n, prefixPi)];
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
    int numberOfNodesWithUnstableQueue = 0;
    for (int i = 0; i < n; ++i) {
        if (oldQ[i] != newQ[i]) {
            ++numberOfNodesWithUnstableQueue;
        }
    }
    return numberOfNodesWithUnstableQueue < 0.1*n;
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
    bool completed = false; // completed when converged and sampled
    do {
        generateNewPackets(n, oldQ, beta, J);
        copy1d(oldQ, newQ, n);

        transmitPackets(n, prefixP, oldQ, newQ);
        copy1d(newQ, oldQ, n);

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
void initNewMemory2d(int n, int m, T*** A) {
    *A = new T*[n];
    for (int i = 0; i < n; ++i) {
        (*A)[i] = new T[m];
    }
}

template <typename T>
void replaceArrayWithPrefixSum(int n, T* A) {
    for (int i = 1; i < n; ++i) {
        A[i] += A[i - 1];
    }
}

void Lsolver::computePrefixP(int n, const Graph *g, double **prefixP) {
    g->copyTransitionMatrix(prefixP);
    for (int i = 0; i < n; ++i) {
        replaceArrayWithPrefixSum(n, prefixP[i]);
    }
}

double computeDistanceFromStationarity(int n, double **P, double *eta) {
    // RMS of (P^T X eta - eta)
    double sumOfSquareError = 0;
    for (int i = 0; i < n; ++i) {
        double error = -eta[i];
        for (int j = 0; j < n; ++j) {
            error += (P[j][i] * eta[j]);
        }
        sumOfSquareError += error * error;
    }
    return sqrt(sumOfSquareError/n);
}

double Lsolver::computeStationaryState(
        const Graph *g, const double *b, double **eta) {

    int n = g->getNumVertex();

    auto T_samp = 4*log(n) / (k*k*e2*e2);

    double *J = new double[n];
    computeJ(n, b, J);

    double **prefixP = NULL;
    initNewMemory2d(n, n, &prefixP);

    // prefixP is replacing every row of P with prefixSum of that row
    // with this pickRandomNeighbor to transmit msg can be done in log(n)
    computePrefixP(n, g, prefixP);

    int *cnt  = new int[n];

    // since transmission is concurrent, we need to hold 2 values
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
    auto error = computeDistanceFromStationarity(n, prefixP, *eta);
    del(n, prefixP);

    std::cerr << "Stationarity distance: " << error << std::endl;
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
