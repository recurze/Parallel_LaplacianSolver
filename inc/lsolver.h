#ifndef LSOLVER_H
#define LSOLVER_H

#include "graph.h"
#include "sampler.h"

#include <vector>

class Lsolver {
public:
    Lsolver(const Graph& g) {
        initGraph(g);
    }

    Lsolver(const Graph& g, const std::vector<double>& b) {
        initGraph(g);
        computeJ(b);
    }

    std::vector<double> solve();
    std::vector<double> solve(const std::vector<double>& b) {
        computeJ(b);
        return solve();
    }
private:
    int n;

    std::vector<double> d;
    std::vector<double> J;
    std::vector<double> eta;

    double beta;
    double b_sink;

    std::vector<int> Q;
    std::vector<int> cnt;

    std::vector<Sampler> sampler;

    void initGraph(const Graph& g);
    void computeJ(const std::vector<double>& b);

    void computeStationarityState();

    void step();
    void serial();
    void parallel();
    void estimateEta();

    std::vector<double> computeX();
};

#endif
