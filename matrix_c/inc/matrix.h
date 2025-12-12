#ifndef MATRIX_LIBRARY_H
#define MATRIX_LIBRARY_H
#include "stdint.h"


typedef struct {
	uint32_t m;
	uint32_t n;
	double det;
	double* data;
} matrix;



matrix* matrix_new(uint32_t m, uint32_t n);
void matrix_del(matrix* mat);


#endif //MATRIX_LIBRARY_H