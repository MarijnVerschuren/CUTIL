#include "matrix.h"

#include <stdlib.h>


matrix* matrix_new(uint32_t m, uint32_t n) {
	matrix* mat = (matrix*)malloc(sizeof(matrix));
	mat->m = m;
	mat->n = n;
	mat->data = (double*)malloc(sizeof(double) * m * n);
	return mat;
}

void matrix_del(matrix* mat) {
	free(mat->data);
	free(mat);
}
