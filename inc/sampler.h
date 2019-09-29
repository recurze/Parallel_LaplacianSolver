#ifndef SAMPLER_H
#define SAMPLER_H

#include "rng.h"

#include <vector>

// Source: http://www.keithschwarz.com/darts-dice-coins/
class Sampler {
public:
    Sampler(const std::vector< std::pair<int, double> >& p);

    int generate() {
        int col = (int) (random_double*n);
        int ind = (random_double < prob[col]) ? col : alias[col];
        return objects[ind];
    }
private:
    int n;
    std::vector<int> objects;

    std::vector<double> prob;
    std::vector<double> alias;

    void init(std::vector<double>& p);
};

#endif
