
#include "mathutil.h"
#include <stdint.h>

double fRand(double a, double b) {
	double f = (double)rand() / RAND_MAX;
	return a + f * (b - a);
}

inline int fast_rand() {
	static uint32_t seed = 0;
	seed = (214013 * seed + 2531011);
	return (seed >> 16) & 0x7fff;
}

#define MY_RANDOM_MAX 0x7fff

double rand_uniform() {
	return (double)fast_rand() / MY_RANDOM_MAX;
}
