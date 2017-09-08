#pragma once

/* Vector operations */

#include <functional>
#include <algorithm>
#include <vector>

// dst = src1 + src2
template<typename T>
inline std::vector<T>& vec_add(const std::vector<T>& src1, const std::vector<T>& src2, std::vector<T>& dst) {
    std::transform(src1.begin(), src1.end(), src2.begin(), dst.begin(), std::plus<T>());
}

// dst = src1 * src2
template<typename T>
inline std::vector<T>& vec_mul(const std::vector<T>& src1, const std::vector<T>& src2, std::vector<T>& dst) {
    std::transform(src1.begin(), src1.end(), src2.begin(), dst.begin(), std::multiplies<T>());
}

// dst = src * scalar
template<typename T>
inline std::vector<T>& vec_mul(const std::vector<T>& src, std::vector<T>& dst, T scalar) {
    std::transform(src.begin(), src.end(), dst.begin(), std::bind1st(std::multiplies<T>(), scalar));
    return dst;
}

// dst = src / scalar
template<typename T>
inline std::vector<T>& vec_div(const std::vector<T>& src, std::vector<T>& dst, T scalar) {
    std::transform(src.begin(), src.end(), dst.begin(), std::bind2nd(std::divides<T>(), scalar));
    return dst;
}

template<typename T>
inline T vec_max(const std::vector<T>& v) {
    return *std::max_element(v.begin(), v.end());
}

inline double vec_sum(const std::vector<double>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0);
}

inline void vec_reset(std::vector<double>& v) {
    std::fill(v.begin(), v.end(), 0.0);
}

inline void vec_reset(std::vector<int>& v) {
    std::fill(v.begin(), v.end(), 0);
}