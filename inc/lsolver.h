#ifndef LSOLVER_H
#define LSOLVER_H

#include "graph.h"

class Lsolver {
public:
    Lsolver(double _e1 = 0.05, double _e2 = 0.05, double _k = 0.1):
        k(_k), e1(_e1), e2(_e2) { }

    void solve(const Graph *g, const double *b, double **x);
private:
    double k;
    double e1;
    double e2;

    void computeJ(int n, const double *b, double *J);
    void computePrefixP(int n, const Graph *g, double **prefixP);
    double computeZstar(int n, const double *eta, const double *d);

    double computeQueueOccupancyProbabilityAtStationarity(
            const Graph *g, const double *b, double **eta);

    void updateCnt(int n, int *Q, int *cnt);
    bool hasConverged(int n, int *inQ);

    void generateNewPackets(int n, int *Q, double beta, const double *J);

    void transmitPackets(int n, double **prefixP, int *oldQ, int *newQ);
    void transmitToRandomNeighbor(int n, double *prefixPi, int *Q, int qid);

    void estimateQueueOccupancyProbability(
            int n, double **P, int *cnt, int *Q, int *inQ,
            double beta, const double *J, double T_samp, double *eta);

    void computeCanonicalSolution(
            const Graph *g, const double *b,
            double *eta, double beta, double **x);
};

#endif
