#include "lsolver.h"

#include "omp.h"
#include <cmath>
#include <random>
#include <cassert>
#include <numeric>
#include <iostream>
#include <algorithm>

template <typename T>
inline T max(int n, const T *a) {
    return *std::max_element(a, a + n);
}

template <typename T>
inline T sum(int n, const T *a) {
    return std::accumulate(a, a + n, (T) 0);
}

template <typename T>
void addArray(int n, T *a, const T *b) {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        a[i] += b[i];
    }
}

template <typename T>
void fill(int n, T *a, const T& x) {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        a[i] = x;
    }
}

template <typename T>
void initNewMemory2d(int n, T*** A) {
    *A = new T*[n];
    for (int i = 0; i < n; ++i) {
        (*A)[i] = new T[n];
    }
}

template <typename T>
void del(int n, T **P) {
    for (int i = 0; i < n; ++i) {
        delete[] P[i];
    }
    delete[] P;
}

void Lsolver::solve(const Graph *g, const double *b, double **x) {
    auto start_time = omp_get_wtime();

    double *eta = NULL;
    auto beta = computeEtaAtStationarity(g, b, &eta);

    assert(eta != NULL);
    assert(beta > 0);

    computeCanonicalSolution(g, b, eta, beta, x);
    delete[] eta; eta = NULL;

    auto end_time = omp_get_wtime();

    std::cerr << "Beta: " << beta
              << "\nTime: " << end_time - start_time
              << "\n";
}

void Lsolver::computeJ(int n, const double *b, double *J) {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        J[i] = -b[i]/b[n - 1];
    }
}

static unsigned long x = 123456789;
static unsigned long y = 362436069;
static unsigned long z = 521288629;

unsigned long xorshf96() {
    unsigned long t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;

    return z;
}

const unsigned long long int MAX = 18446744073709551615ULL;
inline bool trueWithProbability(double p) {
    auto rand = (double) xorshf96()/MAX;
    return rand <= p;
}

void Lsolver::generateNewPackets(
        int n, int *Q, double beta, const double *J) {
#pragma omp parallel for
    for (int i = 0; i < n - 1; ++i) {
        if (trueWithProbability(beta * J[i])) {
            ++Q[i];
        }
    }
}

template <typename T>
int upper_bound(const T *a, int n, const T& x) {
    int l = 0;
    int h = n;
    while (l < h) {
        int m = l + (h - l)/2;
        if (x >= a[m]) {
            l = m + 1;
        } else {
            h = m;
        }
    }
    return l;
}

// Rand lands in [prefixPi[i - 1], prefixPi[i]) with probability P[i]
int Lsolver::pickRandomNeighbor(int n, const double *prefixPi) {
    auto rand = (double) xorshf96()/MAX;
    return upper_bound(prefixPi, n, rand);
}

void Lsolver::transmitPackets(
        int n, double **prefixP, int *Q, int *inQ) {
#pragma omp parallel for
    for (int i = 0; i < n - 1; ++i) {
        if (Q[i] > 0) {
            --Q[i];
#pragma omp atomic
            ++inQ[pickRandomNeighbor(n, prefixP[i])];
        }
    }
}

void Lsolver::updateCnt(int n, const int *Q, int *cnt) {
#pragma omp parallel for
    for (int i = 0; i < n - 1; ++i) {
        if (Q[i] > 0) {
            ++cnt[i];
        }
    }
}

void Lsolver::estimateEta(
        int n, double **prefixP, int *cnt, int *Q, int *inQ,
        double beta, const double *J, double T, double *eta) {

    fill(n, Q, 0);
    fill(n, cnt, 0);

    for (int t = 0; t < T; ++t) {
        fill(n, inQ, 0);

        generateNewPackets(n, Q, beta, J);
        transmitPackets(n, prefixP, Q, inQ);
        addArray(n, Q, inQ);

        updateCnt(n, Q, cnt);
    }

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        eta[i] = cnt[i]/T;
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
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        replaceArrayWithPrefixSum(n, prefixP[i]);
    }
}

double Lsolver::computeEtaAtStationarity(
        const Graph *g, const double *b, double **eta) {

    rng.seed(std::random_device{}());

    int n = g->getNumVertex();
    auto T_samp = 64*3*n + 4*log(n)/(k*e2);

    double *J = new double[n];
    computeJ(n, b, J);

    double **prefixP = NULL;
    initNewMemory2d(n, &prefixP);

    // prefixP is each row of P replaced with prefixSum array of that row
    // to pickNeighbor in O(logN)
    computePrefixP(n, g, prefixP);

    int *cnt  = new int[n];

    // since transmission is concurrent, we need separate inbox
    int *Q = new int[n];
    int *inQ = new int[n];

    *eta = new double[n];
    double max_eta = 0;

    // choose beta such that it's less that beta* (which we don't know)
    // but not too small else packets won't be generated
    // So start with INF or at least
    double beta = 10;
    do {
        beta /= 2;

        estimateEta(n, prefixP, cnt, Q, inQ, beta, J, T_samp, *eta);

        max_eta = max(n, *eta);
    } while (max_eta > 0.75*(1 - e1 - e2) and beta > 0);

    del(n, prefixP);
    delete[] J; J = NULL;
    delete[] Q; Q = NULL;
    delete[] cnt; cnt = NULL;
    delete[] inQ; inQ = NULL;

    return beta;
}

void Lsolver::computeCanonicalSolution(
        const Graph *g, const double *b,
        double *eta, double beta, double **x) {

    int n = g->getNumVertex();

    double* d = new double[n];
    g->copyDegreeMatrix(d);

    *x = new double[n];

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        (*x)[i] = (-b[n - 1]/beta) * (eta[i]/d[i]);
    }

    // centering for canonical solution
    auto avg_x = sum(n, *x)/n;
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        (*x)[i] -= avg_x;
    }

    delete[] d;
}
