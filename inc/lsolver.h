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
    void computeAliasAndProb(
            int n, const Graph *g, double **alias, double **prob);

    double computeEtaAtStationarity(
            const Graph *g, const double *b, double **eta);

    void updateCnt(int n, const int *Q, int *cnt);

    void generateNewPackets(int n, int *Q, double beta, const double *J);

    int pickRandomNeighbor(int n, const double *alias, const double *prob);
    void transmitPackets(
            int n, double **alias, double **prob, int *oldQ, int *newQ);

    void estimateEta(
            int n, double **alias, double **prob, int *cnt, int *Q,
            int *inQ, double beta, const double *J, double *eta);

    void computeCanonicalSolution(
            const Graph *g, const double *b,
            double *eta, double beta, double **x);
};

#endif
