#include "matrix.hpp"
#include "matrix_int.hpp"


namespace MAT {
/*!<
 * Property classes
 * */
f64_t& det_t::operator =(const f64_t val) {
	return this->_det = val;
}


/*!<
 * matrix
 * */
matrix::matrix(const uint32_t n, const uint32_t m) {
	this->n = n;
	this->m = m;
	this->data = (f64_t*)malloc(sizeof(f64_t) * n * m);
	// TODO: aligned_alloc
}

matrix::~matrix(void) {
	free(this->data);
}

void matrix::print(void) const {
	for (uint32_t i = 0; i < this->n; i++) {
		printf("[");
		for (uint32_t j = 0; j < this->m; j++) {
			printf("%5.2f ", this->data[i * this->m + j]);
		}
		printf("]\n");
	}
	printf("\n");
}

matrix& matrix::reduce(void) {
	//matrix_swap_rows(this->data, this->m, 0, this->n - 1);
	matrix_add_rows(this->data, this->m, 2, 0, -2.0f);
	return *this;
}

double& matrix::operator()(const uint32_t i, const uint32_t j) {
	return this->data[i * m + j];
}




/*!<
 * square matrix
 * */
sq_matrix& sq_matrix::L_matrix(void) {
	if (this->LU[0] == nullptr) { LU_decomposition(); }
	return *this->LU[0];
}

sq_matrix& sq_matrix::U_matrix(void) {
	if (this->LU[1] == nullptr) { LU_decomposition(); }
	return *this->LU[1];
}


void sq_matrix::LU_decomposition(void) {
	// TODO
}

void sq_matrix::new_RRWS(void) {
	// TODO: CHANGE!!!
	RRWS = new matrix(this->n, this->m * 2);
	for (uint32_t i = 0; i < n; i++) {
		memcpy(&RRWS->data[i*m*2], &this->data[i * m], sizeof(f64_t) * m);
		RRWS->data[i*m*2 + m + i] = 1.0;	// identity matrix
	}
}

}