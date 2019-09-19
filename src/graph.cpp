#include "graph.h"

#include <numeric>

template <typename T>
inline T sum(const std::vector<T>& a) {
    return std::accumulate(a.begin(), a.end(), (T) 0);
}

Graph::Graph(const std::vector< std::vector<double> >& A) {
    n = (int) A.size();
    D.resize(n);
    P.resize(n);
    for (int i = 0; i < n; ++i) {
        D[i] = sum(A[i]);

        P[i].resize(n);
        for (int j = 0; j < n; ++j) {
            P[i][j] = A[i][j]/D[i];
        }
    }
}
