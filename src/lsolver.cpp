#include "rng.h"
#include "lsolver.h"

#include <omp.h>

#include <cmath>
#include <chrono>
#include <cassert>
#include <numeric>
#include <iostream>
#include <algorithm>

template <typename T>
inline T max(const std::vector<T>& a) {
    return *std::max_element(a.begin(), a.end());
}

template <typename T>
inline T sum(const std::vector<T>& a) {
    return std::accumulate(a.begin(), a.end(), (T) 0);
}

void Lsolver::initGraph(const Graph& g) {
    n = g.getNumVertex();
    d = g.getDegreeMatrix();

    sampler.reserve(n);
    for (int i = 0; i < n; ++i) {
        auto v(g.getNeighbors(i));
        for (auto& j: v) {
            j.second /= d[i];
        }
        sampler.push_back(Sampler(v));
    }
}

void Lsolver::computeJ(const std::vector<double>& b) {
    b_sink = b.back();
    J.resize(n);
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        J[i] = -b[i]/b_sink;
    }
}

std::vector<double> Lsolver::solve() {
    auto start = std::chrono::steady_clock::now();

    computeStationarityState();
    auto x = computeX();

    auto finish = std::chrono::steady_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<
        std::chrono::duration<double> >(finish - start).count();
    std::cerr << "Solve Time: " << elapsed_seconds << '\n';

    return x;
}

inline bool trueWithProbability(double p) {
    return random_double <= p;
}

void Lsolver::step() {
    parallel();
    //serial();
}

#ifndef N_THREADS
#define N_THREADS 16
#endif

void Lsolver::parallel() {
    std::vector<int> inQ[N_THREADS];

#pragma omp parallel for
    for (int i = 0; i < N_THREADS; ++i) {
        inQ[i].resize(n, 0);
    }

#pragma omp parallel for
    for (int i = 0; i < n - 1; ++i) {
        Q[i] += trueWithProbability(beta * J[i]);
        for (int k = 3; Q[i] and k; --k) {
            --Q[i];
            ++cnt[i];
            ++inQ[omp_get_thread_num()][sampler[i].generate()];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        int x = 0;
        for (int p = 0; p < N_THREADS; ++p) {
            x += inQ[p][i];
        }
        Q[i] += x;
    }
}

void Lsolver::serial() {
    std::vector<int> inQ(n, 0);
    for (int i = 0; i < n - 1; ++i) {
        Q[i] += trueWithProbability(beta * J[i]);
        for (int k = 3; Q[i] and k; --k) {
            --Q[i];
            ++cnt[i];
            ++inQ[sampler[i].generate()];
        }
    }
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        Q[i] += inQ[i];
    }
}

const double EPS = 1e-3;
const int MAX_EPOCHS = 1000;
const int LENGTH_OF_EPOCH = 1500;

void Lsolver::estimateEta() {
    int epoch = 0;
    double oldC = 0, newC = 0;
    do {
        ++epoch;
        oldC = newC;
        for (int t = 0; t < LENGTH_OF_EPOCH; ++t) {
            step();
        }
        newC = (double) Q[n - 1]/(1 + sum(Q));
        //std::cerr << newC << '\n';
    } while (fabs(oldC - newC) > EPS and epoch < MAX_EPOCHS);

    int T = epoch * LENGTH_OF_EPOCH;
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        eta[i] = (double) cnt[i]/T;

        Q[i] = 0;
        cnt[i] = 0;
    }
}

void Lsolver::computeStationarityState() {
    omp_set_num_threads(N_THREADS);
    // Can start with any big value, but beta < beta* is below 1
    beta = 1.28;
    eta.resize(n);

    Q.resize(n, 0);
    cnt.resize(n, 0);
    do {
        beta /= 2;
        estimateEta();
    } while (max(eta) > 0.75 and beta > 0);
}

std::vector<double> Lsolver::computeX() {
    std::vector<double> x(n);
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        x[i] = (-b_sink/beta) * (eta[i]/d[i]);
    }
    return x;

    // centering for canonical solution
    //auto avg_x = sum(x)/n;
    //for (int i = 0; i < n; ++i) {
    //    x[i] -= avg_x;
    //}
}
