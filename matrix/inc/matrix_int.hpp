//
// Created by marijn on 12/8/25.
//
#ifndef MATRIX_INT_H
#define MATRIX_INT_H
#include "types.h"


/*!
 * matrix_row_ops.asm
 */
extern "C" void matrix_swap_rows(f64_t* data, u32_t m, u32_t d, u32_t s);
extern "C" void matrix_add_rows(f64_t* data, u32_t m, u32_t d, u32_t s, f64_t scalar);
extern "C" void matrix_scale_row(f64_t* data, u32_t m, u32_t r, f64_t scalar);


/*!
 * matrix_REF/RREF.asm
 */
extern "C" void matrix_REF_64(f64_t* data, u32_t m, u32_t n);
extern "C" void matrix_RREF_64(f64_t* data, u32_t m, u32_t n);
extern "C" void matrix_REF_32(f32_t* data, u32_t m, u32_t n);
extern "C" void matrix_RREF_32(f32_t* data, u32_t m, u32_t n);


/*!
 * matrix_math.asm
 */
extern "C" void matrix_add_64(f64_t* dst, f64_t* a, f64_t* b, u32_t cnt);
extern "C" void matrix_add_32(f32_t* dst, f32_t* a, f32_t* b, u32_t cnt);
extern "C" void matrix_sub_64(f64_t* dst, f64_t* a, f64_t* b, u32_t cnt);
extern "C" void matrix_sub_32(f32_t* dst, f32_t* a, f32_t* b, u32_t cnt);
extern "C" void matrix_scale_64(f64_t* dst, f64_t* src, f64_t scalar, u32_t cnt);
extern "C" void matrix_scale_32(f32_t* dst, f32_t* src, f32_t scalar, u32_t cnt);


/*!
 * matrix_math.cpp
 */
void matrix_mul_64(f64_t* C, f64_t* A, f64_t* B, u32_t m, u32_t n, u32_t p);	// C = A*B (A[mXn], b[nXp])
void matrix_mul_32(f32_t* C, const f32_t* A, f32_t* B, u32_t m, u32_t n, u32_t p);	// C = A*B (A[mXn], b[nXp])


//void matrix_REF_C(double* data, uint32_t n, uint32_t m) {
//	f64_t pivot;
//	for (uint32_t j, i = 0; i < n; i++) {
//		for (j = i; j < n; j++) {
//			pivot = data[j * m + i];
//			if (pivot != 0.0f) { break; }
//		}
//		if (pivot == 0.0f) { continue; }
//		if (i != j) { matrix_swap_rows(data, m, i, j); }
//		for (j = i + 1; j < n; j++) {
//			matrix_add_rows(data, m, j, i, -data[j * m + i]/pivot);
//		}
//	}
//}
//
//void matrix_RREF_C(double* data, uint32_t n, uint32_t m) {
//	f64_t pivot;
//	for (uint32_t j, i = 0; i < n; i++) {
//		for (j = i; j < n; j++) {
//			pivot = data[j * m + i];
//			if (pivot != 0.0f) { break; }
//		}
//		if (pivot == 0.0f) { continue; }
//		if (i != j) { matrix_swap_rows(data, m, i, j); }
//
//		matrix_scale_row(data, m, i, 1/pivot);
//		pivot = 1.0f;
//
//		for (j = 0; j < n; j++) {
//			if (j == i) { continue; }
//			matrix_add_rows(data, m, j, i, -data[j * m + i]/pivot);
//		}
//	}
//}

// TODO: PLU decomposition
// NOTE: square only
//extern "C" void matrix_REF_LU(double* data, uint32_t dim, double* L);

#endif //MATRIX_INT_H
