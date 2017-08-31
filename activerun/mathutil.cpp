
#include "mathutil.h"
#include <stdint.h>

uint32_t seed = 0;

inline int fast_rand() {
	seed = (214013 * seed + 2531011);
	return (seed >> 16) & 0x7fff;
}

#ifndef CUSTOM_RANDOM

#define MY_RANDOM_MAX RAND_MAX

double rand_uniform() {
	return (double)rand() / MY_RANDOM_MAX;
}

void set_random_seed(int rseed) {
    srand(rseed);
}

#else

#define MY_RANDOM_MAX 0x7fff

double rand_uniform() {
    return (double)fast_rand() / MY_RANDOM_MAX;
}

void set_random_seed(int rseed) {
    seed = rseed;
}


#endif // !CUSTOM_RANDOM

