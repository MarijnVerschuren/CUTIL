//
// Created by marijn on 12/8/25.
//
#ifndef MATRIX_INT_H
#define MATRIX_INT_H
#include "types.h"

/*!
 * matrix_row_ops.asm
 */
extern "C" void matrix_swap_rows(f64_t* data, uint32_t m, uint32_t d, uint32_t s);
extern "C" void matrix_add_rows(f64_t* data, uint32_t m, uint32_t d, uint32_t s, f64_t scalar);
extern "C" void matrix_scale_row(f64_t* data, uint32_t m, uint32_t r, f64_t scalar);


//extern "C" void matrix_copy_slice(f64_t* dst, const f64_t* src, uint32_t n, uint32_t m, uint32_t m_start, uint32_t m_end);

// matrix reduce echelon form
extern "C" void matrix_REF_64(f64_t* data, uint32_t n, uint32_t m);
extern "C" void matrix_RREF_64(f64_t* data, uint32_t n, uint32_t m);
extern "C" void matrix_REF_32(f32_t* data, uint32_t n, uint32_t m);
extern "C" void matrix_RREF_32(f32_t* data, uint32_t n, uint32_t m);
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
