#ifndef LSOLVER_H
#define LSOLVER_H

#include "graph.h"

class Lsolver {
public:
    Lsolver(double _e1 = 0.1, double _e2 = 0.1, double _k = 0.1):
        k(_k), e1(_e1), e2(_e2) { }

    void solve(const Graph *g, const double *b, double *x);
private:
    double k;
    double e1;
    double e2;

    void computeJ(int n, const double *b, double *J);
    double computeZstar(int n, const double *eta, const double *d);

    void estimateQueueOccupancyProbability(
            int n, double **P, double beta,
            const double *J, double T_samp, double *eta);
    double computeStationaryState(
            const Graph *g, const double *b, double *eta);

    void computeCanonicalSolution(
            const Graph *g, const double *b,
            double *eta, double beta, double *x);
};

#endif
