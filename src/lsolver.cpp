#include "rng.h"
#include "lsolver.h"

#include <omp.h>

#include <queue>
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

void Lsolver::p1() {
#pragma omp parallel
    {
        auto id = omp_get_thread_num();
        partition[id].reserve(n/N_THREADS);

#pragma omp for
        for (int i = 0; i < n - 1; ++i) {
            partition[id].push_back(i);
        }
    }
}

void Lsolver::p2() {
    std::vector<bool> visited(n, false);
    std::vector<int> a; a.reserve(n - 1);
    auto bfs = [&]() {
        std::queue<int> q; q.push(0);
        visited[0] = true;
        while (not q.empty()) {
            int u = q.front(); q.pop();
            if (u != n - 1) a.push_back(u);
            for (auto [v, w]: adj[u]) {
                if (not visited[v]) {
                    q.push(v);
                    visited[v] = true;
                }
            }
        }
    }; bfs();

#pragma omp parallel
    {
        auto id = omp_get_thread_num();
        partition[id].reserve(n/N_THREADS);

#pragma omp for
        for (int i = 0; i < n - 1; ++i) {
            partition[id].push_back(a[i]);
        }
    }
}

void Lsolver::partitionGraph() {
    p2();
}

void Lsolver::initGraph(const Graph& g) {
    n = g.getNumVertex();
    d = g.getDegreeMatrix();

    adj.resize(n);
    sampler.reserve(n);
    for (int i = 0; i < n; ++i) {
        adj[i] = g.getNeighbors(i);
        for (auto& j: adj[i]) {
            j.second /= d[i];
        }
        sampler.push_back(Sampler(adj[i]));
    }

    partitionGraph();
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

#ifndef N_THREADS
#define N_THREADS 16
#endif

#define LENGTH_OF_EPOCH std::max(n, 2000)
void Lsolver::parallel() {
    std::vector< std::vector<int> > inQ(N_THREADS, std::vector<int>(n, 0));
#pragma omp parallel
    {
        auto id = omp_get_thread_num();
        for (int t = 0; t < LENGTH_OF_EPOCH; ++t) {

#pragma omp for
            for (int i = 0; i < n - 1; ++i) {
                Q[i] += trueWithProbability(beta * J[i]);
                for (int k = 3; Q[i] > 0 and k; --k) {
                    --Q[i];
                    ++cnt[i];
                    ++inQ[id][sampler[i].generate()];
                }
            }

#pragma omp for
            for (int i = 0; i < n; ++i) {
                int x = 0;
                for (int p = 0; p < N_THREADS; ++p) {
                    x += inQ[p][i];
                    inQ[p][i] = 0;
                }
                Q[i] += x;
            }

        }
    }
}

void Lsolver::parallel_partitioned() {
    std::vector<int> inQ(n, 0);
#pragma omp parallel
    {
        auto id = omp_get_thread_num();
        for (int t = 0; t < LENGTH_OF_EPOCH; ++t) {
#pragma omp barrier
            for (int i: partition[id]) {
                Q[i] += trueWithProbability(beta * J[i]);
                for (int k = 3; Q[i] and k; --k) {
                    --Q[i];
                    ++cnt[i];
#pragma omp atomic
                    ++inQ[sampler[i].generate()];
                }
            }

#pragma omp barrier
            for (int i: partition[id]) {
                Q[i] += inQ[i];
                inQ[i] = 0;
            }

#pragma omp single
            {
                Q[n - 1] += inQ[n - 1];
                inQ[n - 1] = 0;
            }
        }
    }
}

void Lsolver::serial() {
    //std::vector<double> inQ(n, 0);
    //for (int i = 0; i < n - 1; ++i) {
    //    for (const auto& j: adj[i]) {
    //        inQ[j.first] += Q[i] * j.second;
    //    }
    //}
    //for (int i = 0; i < n - 1; ++i) {
    //    Q[i] = inQ[i] + beta*J[i];
    //    cnt[i] += Q[i];
    //}
    //Q[n - 1] += inQ[n - 1];

    std::vector<int> inQ(n, 0);
    for (int t = 0; t < LENGTH_OF_EPOCH; ++t) {
        for (int i = 0; i < n - 1; ++i) {
            Q[i] += trueWithProbability(beta * J[i]);
            for (int k = 3; Q[i] and k; --k) {
                --Q[i];
                ++cnt[i];
                ++inQ[sampler[i].generate()];
            }
        }
        for (int i = 0; i < n; ++i) {
            Q[i] += inQ[i];
            inQ[i] = 0;
        }
    }
}

const double EPS = 1e-3;
const int MIN_EPOCHS = 10;
const int MAX_EPOCHS = 100;

double Lsolver::estimateEta(double lastC = 0) {
    int epoch = 0;
    double oldC = 0, newC = 0;

    while (epoch < MIN_EPOCHS and newC < 0.8) {
        ++epoch;
        parallel();
        newC = (double) Q[n - 1]/(1 + sum(Q));
    }

    do {
        ++epoch;

        oldC = newC;
        parallel();
        newC = (double) Q[n - 1]/(1 + sum(Q));

        std::cerr << "C: " << newC << ' '
                  << "Sunk: " << Q[n - 1]
                  << '\n';
        if (beta > 0.001 and fabs(oldC - newC) < EPS and newC > lastC)
            break;
    } while (epoch < MAX_EPOCHS);
    std::cerr << "Epochs: " << epoch << '\n';

    int T = epoch * LENGTH_OF_EPOCH;
    for (int i = 0; i < n; ++i) {
        eta[i] = (double) cnt[i]/T;

        Q[i] = 0, cnt[i] = 0;
    }
    return newC;
}

void Lsolver::computeStationarityState() {
    omp_set_num_threads(N_THREADS);
    // Can start with any big value, but beta < beta* is below 1
    beta = 0.05;
    eta.resize(n);

    double C = 0;

    Q.resize(n, 0);
    cnt.resize(n, 0);
    do {
        beta /= 2;
        std::cerr << "Beta: " << beta << '\n';
        C = estimateEta(C);
    } while ((C < 0.80 or max(eta) > 0.70) and beta >= 0.0001);
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
