#pragma once

#ifndef ARRAYUTIL_H
#define ARRAYUTIL_H

#include "vec.h"
#include <vector>
//#include <immintrin.h>

inline void array_add_double(const double* src, double* dst, unsigned int size) {
	
/*	unsigned int i, j;
	for (i = 0; i < size / 32; i++) {
		_mm256_store_pd(dst + 32 * i + 0, _mm256_addsub_pd(_mm256_load_pd(src + 32 * i + 0), _mm256_load_pd(dst + 32 * i + 0)));
		_mm256_store_pd(dst + 32 * i + 4, _mm256_addsub_pd(_mm256_load_pd(src + 32 * i + 4), _mm256_load_pd(dst + 32 * i + 4)));
		_mm256_store_pd(dst + 32 * i + 8, _mm256_addsub_pd(_mm256_load_pd(src + 32 * i + 8), _mm256_load_pd(dst + 32 * i + 8)));
		_mm256_store_pd(dst + 32 * i +12, _mm256_addsub_pd(_mm256_load_pd(src + 32 * i +12), _mm256_load_pd(dst + 32 * i +12)));
		_mm256_store_pd(dst + 32 * i +16, _mm256_addsub_pd(_mm256_load_pd(src + 32 * i +16), _mm256_load_pd(dst + 32 * i +16)));
		_mm256_store_pd(dst + 32 * i +20, _mm256_addsub_pd(_mm256_load_pd(src + 32 * i +20), _mm256_load_pd(dst + 32 * i +20)));
		_mm256_store_pd(dst + 32 * i +24, _mm256_addsub_pd(_mm256_load_pd(src + 32 * i +24), _mm256_load_pd(dst + 32 * i +24)));
		_mm256_store_pd(dst + 32 * i +28, _mm256_addsub_pd(_mm256_load_pd(src + 32 * i +28), _mm256_load_pd(dst + 32 * i +28)));
	}

	for (j = i * 32; j < size - 4; j += 4) {
		_mm256_store_pd(dst + j, _mm256_addsub_pd(_mm256_load_pd(src + j), _mm256_load_pd(dst + j)));
	}*/

	for (int j = 0; j < size; j++) {
		dst[j] += src[j];
	}
}
/*
inline void array_mul_double(const double* src1, const double* src2, double* dst, unsigned int size) {

	unsigned int i, j;
	for (i = 0; i < size / 32; i++) {
		_mm256_store_pd(dst + 32 * i + 0, _mm256_mul_pd(_mm256_load_pd(src1 + 32 * i + 0), _mm256_load_pd(src2 + 32 * i + 0)));
		_mm256_store_pd(dst + 32 * i + 4, _mm256_mul_pd(_mm256_load_pd(src1 + 32 * i + 4), _mm256_load_pd(src2 + 32 * i + 4)));
		_mm256_store_pd(dst + 32 * i + 8, _mm256_mul_pd(_mm256_load_pd(src1 + 32 * i + 8), _mm256_load_pd(src2 + 32 * i + 8)));
		_mm256_store_pd(dst + 32 * i + 12, _mm256_mul_pd(_mm256_load_pd(src1 + 32 * i + 12), _mm256_load_pd(src2 + 32 * i + 12)));
		_mm256_store_pd(dst + 32 * i + 16, _mm256_mul_pd(_mm256_load_pd(src1 + 32 * i + 16), _mm256_load_pd(src2 + 32 * i + 16)));
		_mm256_store_pd(dst + 32 * i + 20, _mm256_mul_pd(_mm256_load_pd(src1 + 32 * i + 20), _mm256_load_pd(src2 + 32 * i + 20)));
		_mm256_store_pd(dst + 32 * i + 24, _mm256_mul_pd(_mm256_load_pd(src1 + 32 * i + 24), _mm256_load_pd(src2 + 32 * i + 24)));
		_mm256_store_pd(dst + 32 * i + 28, _mm256_mul_pd(_mm256_load_pd(src1 + 32 * i + 28), _mm256_load_pd(src2 + 32 * i + 28)));
	}

	for (j = i * 32; j < size - 4; j += 4) {
		_mm256_store_pd(dst + j, _mm256_mul_pd(_mm256_load_pd(src1 + j), _mm256_load_pd(src2 + j)));
	}

	for (; j < size; j++) {
		dst[j] += src1[j];
	}
}*/


// speed like a wind
inline void array_add(const std::vector<Vec2>& src, std::vector<Vec2>& dst) {
	array_add_double((double*)&src[0], (double*)&dst[0], 2 * src.size());
}

inline void array_add(const std::vector<Vec3>& src, std::vector<Vec3>& dst) {
    array_add_double((double*)&src[0], (double*)&dst[0], 3 * src.size());
}
/*
inline void array_mul(const std::vector<Vec2>& src1, const std::vector<Vec2>& src2, std::vector<Vec2>& dst) {
	array_mul_double((double*)&src1[0], (double*)&src2[0], (double*)&dst[0], 2 * src1.size());
}*/

#endif // !ARRAYUTIL_H

