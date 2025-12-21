//
// Created by marijn on 12/20/25.
//
#include "matrix_int.hpp"

#include <immintrin.h>



// block sizes multiplication TODO: tune
#define BM 64
#define BN 64
#define BK 64


// TODO:
// extern "C" void matrix_add_64(f64_t* dst, f64_t* a, f64_t* b, u32_t cnt);
// extern "C" void matrix_add_32(f32_t* dst, f32_t* a, f32_t* b, u32_t cnt);
// extern "C" void matrix_sub_64(f64_t* dst, f64_t* a, f64_t* b, u32_t cnt);
// extern "C" void matrix_sub_32(f32_t* dst, f32_t* a, f32_t* b, u32_t cnt);
// extern "C" void matrix_scale_64(f64_t* dst, f64_t* src, f64_t scalar, u32_t cnt);
// extern "C" void matrix_scale_32(f32_t* dst, f32_t* src, f32_t scalar, u32_t cnt);


// TODO: inline into method???
void matrix_mul_64(f64_t* C, const f64_t* A, const f64_t* B, u32_t m, u32_t n, u32_t p) {
    u32_t i, j, k, ii, jj, kk;

    __m256d zero = _mm256_setzero_pd();
    for (i = 0; i + 3 < n*p; i += 4) {
        _mm256_storeu_pd(C + i, zero);
    } for (; i < n*p; ++i) {
        C[i] = 0.0;
    }

    for (ii = 0; ii < m; ii += BM) {
        for (kk = 0; kk < n; kk += BK) {
            for (jj = 0; jj < p; jj += BN) {
                u32_t i_max = (ii + BM < m ? ii + BM : m);
                u32_t k_max = (kk + BK < n ? kk + BK : n);
                u32_t j_max = (jj + BN < p ? jj + BN : p);

                for (i = ii; i < i_max; ++i) {
                    for (k = kk; k < k_max; ++k) {
                        f64_t a_val = A[i * n + k];

                        // Pointers for inner-most loop
                        const f64_t* b_ptr = B + k * p + jj;
                        f64_t* c_ptr = C + i * p + jj;

                        u32_t j_inner = 0;

                        // Use AVX2: 4 doubles per iteration
                        for (; j_inner + 3 < j_max - jj; j_inner += 4) {
                            __m256d c_vec = _mm256_loadu_pd(c_ptr + j_inner);
                            __m256d b_vec = _mm256_loadu_pd(b_ptr + j_inner);
                            __m256d a_vec = _mm256_set1_pd(a_val); // broadcast a_val
                            c_vec = _mm256_fmadd_pd(a_vec, b_vec, c_vec); // c += a*b
                            _mm256_storeu_pd(c_ptr + j_inner, c_vec);
                        }

                        // handle tail
                        for (; j_inner < j_max - jj; ++j_inner) {
                            c_ptr[j_inner] += a_val * b_ptr[j_inner];
                        }
                    }
                }
            }
        }
    }
}

void matrix_mul_32(f32_t* C, const f32_t* A, const f32_t* B, u32_t m, u32_t n, u32_t p) {
    u32_t i, j, k, ii, jj, kk;

    const __m256 zero = _mm256_setzero_ps();
    for (i = 0; i + 7 < n*p; i += 8) {
        _mm256_storeu_ps(C + i, zero);
    } for (; i < n*p; ++i) {
        C[i] = 0.0;
    }

    for (ii = 0; ii < m; ii += BM) {
        for (kk = 0; kk < n; kk += BK) {
            for (jj = 0; jj < p; jj += BN) {
                u32_t i_max = (ii + BM < m ? ii + BM : m);
                u32_t k_max = (kk + BK < n ? kk + BK : n);
                u32_t j_max = (jj + BN < p ? jj + BN : p);

                for (i = ii; i < i_max; ++i) {
                    for (k = kk; k < k_max; ++k) {
                        f32_t a_val = A[i * n + k];

                        // Pointers for inner-most loop
                        const f32_t* b_ptr = B + k * p + jj;
                        f32_t* c_ptr = C + i * p + jj;

                        u32_t j_inner = 0;

                        // Use AVX2: 4 doubles per iteration
                        for (; j_inner + 7 < j_max - jj; j_inner += 8) {
                            __m256 c_vec = _mm256_loadu_ps(c_ptr + j_inner);
                            __m256 b_vec = _mm256_loadu_ps(b_ptr + j_inner);
                            __m256 a_vec = _mm256_set1_ps(a_val); // broadcast a_val
                            c_vec = _mm256_fmadd_ps(a_vec, b_vec, c_vec); // c += a*b
                            _mm256_storeu_ps(c_ptr + j_inner, c_vec);
                        }

                        // handle tail
                        for (; j_inner < j_max - jj; ++j_inner) {
                            c_ptr[j_inner] += a_val * b_ptr[j_inner];
                        }
                    }
                }
            }
        }
    }
}