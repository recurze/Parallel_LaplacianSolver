#ifndef RNG_H
#define RNG_H

class Rng {
public:
    Rng();

    double random_double();
private:
    unsigned long x = 123456789;
    unsigned long y = 362436069;
    unsigned long z = 521288629;

    unsigned long xorshf96();
};
#endif
