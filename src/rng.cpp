#include "rng.h"

unsigned long x = 123456789;
unsigned long y = 362436069;
unsigned long _z = 521288629;
inline unsigned long xorshf96() {
    unsigned long t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

    t = x;
    x = y;
    y = _z;
    _z = t ^ x ^ y;

    return _z;
}

// Source: http://www.cse.yorku.ca/~oz/marsaglia-rng.html
#define znew (z=36969*(z&65535)+(z>>16))
#define wnew (w=18000*(w&65535)+(w>>16))
#define MWC  ((znew<<16)+wnew)
#define CONG (jcong=69069*jcong+1234567)
#define SHR3 (jsr^=(jsr<<17), jsr^=(jsr>>13), jsr^=(jsr<<5))
#define KISS ((MWC^CONG)+SHR3)

static unsigned long z = 12345;
static unsigned long w = 65435;
static unsigned long jsr = 34221;
static unsigned long jcong = 12345;
const unsigned long MAX = 18446744073709551615ULL;

double Rng::random_double() {
    return (double) KISS/MAX;
}

