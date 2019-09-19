#ifndef SAMPLER_H
#define SAMPLER_H

#include "rng.h"

#include <vector>

class Sampler {
public:
    Sampler(std::vector<double>& p);

    int generate();
private:
    int n;
    std::vector<double> prob;
    std::vector<double> alias;

    Rng rng;
    void init(std::vector<double>& p);
};

#endif
