#ifndef RNG_H
#define RNG_H

#define MAX 18446744073709551615ULL

// Source: http://www.cse.yorku.ca/~oz/marsaglia-rng.html
static unsigned long z = 12345;
static unsigned long w = 65435;
static unsigned long jsr = 34221;
static unsigned long jcong = 12345;

#define znew (z=36969*(z&65535)+(z>>16))
#define wnew (w=18000*(w&65535)+(w>>16))
#define MWC  ((znew<<16)+wnew)
#define CONG (jcong=69069*jcong+1234567)
#define SHR3 (jsr^=(jsr<<17), jsr^=(jsr>>13), jsr^=(jsr<<5))
#define KISS ((MWC^CONG)+SHR3)

#define random_double ((double) KISS/MAX)

#endif
