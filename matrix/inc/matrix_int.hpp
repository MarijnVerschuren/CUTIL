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



//extern "C" void matrix_copy_slice(f64_t* dst, const f64_t* src, uint32_t n, uint32_t m, uint32_t m_start, uint32_t m_end);

// matrix reduce echelon form
extern "C" void matrix_REF(f64_t* data, uint32_t n, uint32_t m);
//void matrix_REF(double* data, uint32_t n, uint32_t m) {
//	double pivot, scalar;
//
//	// TODO try to remove
//	double* pivoting_row;
//	uint32_t pivoting_col;
//
//	for (uint32_t i = 0; i < n; i++) {
//		pivoting_row = &data[i * m];
//		pivoting_col = ((i < m) ? i : m);
//		do {
//			pivot = pivoting_row[pivoting_col];
//		} while (pivot == 0 && ++pivoting_col < m);
//		for (uint32_t i2 = i + 1; i2 < m; i2++) {
//			scalar = data[i2 * m + pivoting_col] / pivot;
//			for (uint32_t j = i; j < m; j++) {
//				data[i2 * m + j] -= scalar * pivoting_row[j];
//			}
//		}
//	}
//}

// TODO: PLU decomposition
// NOTE: square only
//extern "C" void matrix_REF_LU(double* data, uint32_t dim, double* L);

#endif //MATRIX_INT_H
