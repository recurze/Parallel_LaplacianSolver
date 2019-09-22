#ifndef SAMPLER_H
#define SAMPLER_H

#include "rng.h"

#include <vector>

// Source: http://www.keithschwarz.com/darts-dice-coins/
class Sampler {
public:
    Sampler(const std::vector<double>& p);

    int generate();
private:
    int n;
    std::vector<double> prob;
    std::vector<double> alias;

    Rng rng;
    void init(const std::vector<double>& p);
};

#endif
