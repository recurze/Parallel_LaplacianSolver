#ifndef LSOLVER_H
#define LSOLVER_H

#include "graph.h"

class Lsolver {
public:
    Lsolver(): e1(0.1), e2(0.1) { }
    Lsolver(double _e1, double _e2): e1(_e1), e2(_e2) { }

    void solve(const Graph *g, const double *b, double *x);
private:
    double e1;
    double e2;
    double computeT_hit(const Graph *g);

    void estimateQueueOccupancyProbability(
            double **P, const double beta, const double *J,
            const double T_mix, const double T_samp, double *eta);
    double computeStationaryState(
            const Graph *g, const double *b, double t_hit, double *eta);

    void computeCanonicalSolution(
            const Graph *g, const double *b, double *eta, double beta, double *x);
};

#endif
