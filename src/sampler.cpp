#include "sampler.h"

Sampler::Sampler(std::vector<double>& p) {
    n = (int) p.size();
    prob.resize(n);
    alias.resize(n);
    init(p);
}

void Sampler::init(std::vector<double>& P) {
    std::vector<int> small; small.reserve(n);
    std::vector<int> large; large.reserve(n);

    for (int i = 0; i < n; ++i) {
        P[i] *= n;
        if (P[i] < 1) {
            small.push_back(i);
        } else {
            large.push_back(i);
        }
    }

    while (not small.empty() and  not large.empty()) {
        auto less = small.back(); small.pop_back();
        auto more = large.back(); large.pop_back();

        prob[less] = P[less];
        alias[less] = more;

        P[more] -= (1 - P[less]);

        if (P[more] < 1) {
            small.push_back(more);
        } else {
            large.push_back(more);
        }
    }

    for (auto &i: large) prob[i] = 1;
    for (auto &i: small) prob[i] = 1;
}

int Sampler::generate() {
    int col = (int) (rng.random_double()*n);
    return (rng.random_double() < prob[col]) ? col : alias[col];
}
