#pragma once

#ifndef ARRAYUTIL_H
#define ARRAYUTIL_H

#include "vec.h"



// speed like a wind
void array_add(const double* src, double* dst, size_t len) {
    const double* end = src + len;
    while (src < end) {
        *(dst++) += *(src++);
    }
}

void array_add(Vec2* src, Vec2* dst, size_t len) {
    array_add((double*)&src[0], (double*)&dst[0], len + len);
}

void array_mul(double* src, double d, size_t len) {
    double* end = src + len;
    while (src < end) {
        *(src++) *= d;
    }
}

void array_mul(Vec2* src, double d, size_t len) {
    array_mul((double*)&src[0], d, len + len);
}

void array_mul(const double* src, double* dst, size_t len) {
    const double* end = src + len;
    while (src < end) {
        *(dst++) = *(dst++) * (*src);
    }
}
void array_mul(const double* src, double* dst, double d, size_t len) {
    const double* end = src + len;
    while (src < end) {
        *(dst++) = *(src++) * d;
    }
}
void array_mul_add(const double* src, double* dst, double d, size_t len) {
    const double* end = src + len;
    while (src < end) {
        *(dst++) += *(src++) * d;
    }
}

void array_div(const double* src, double* dst, double d, size_t len) {
    const double* end = src + len;
    while (src < end) {
        *(dst++) = d / *(src++);
    }
}

template<typename T>
struct FixArray {
    unsigned int size;
    T* data;

    inline void clear() {
        size = 0;
    }

    T operator[](int idx)const {
        return data[idx];
    }
    T& operator[](int idx) {
        return data[idx];
    }

    void push_back(T t) {
        data[size++] = t;
    }

    T* begin() {
        return data;
    }
    T* end() {
        return data + size;
    }
};

#endif // !ARRAYUTIL_H

