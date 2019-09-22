#ifndef LSOLVER_H
#define LSOLVER_H

#include "graph.h"
#include "sampler.h"

#include <vector>

class Lsolver {
public:
    Lsolver(const Graph *g);
    Lsolver(const Graph *g, const std::vector<double>& b);

    void solve(std::vector<double>& x);
    void solve(const std::vector<double>& b, std::vector<double>& x);
private:
    int n;

    std::vector<double> d;
    std::vector<double> J;
    std::vector<double> eta;

    double beta;
    double b_sink;

    std::vector< std::vector<double> > P;

    std::vector<Sampler> sampler;

    const double k = 0.1;
    const double e1 = 0.1;
    const double e2 = 0.1;

    void initGraph(const Graph *g);

    void computeJ(const std::vector<double>& b);

    void computeStationarityState();

    void estimateEta();

    void computeCanonicalSolution(std::vector<double>& x);
};

#endif
