
#include "mathutil.h"

double fRand(double a, double b) {
	double f = (double)rand() / RAND_MAX;
	return a + f * (b - a);
}

double rand_uniform() {
	return (double)rand() / RAND_MAX;
}
