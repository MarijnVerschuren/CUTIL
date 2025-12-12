#include "matrix.hpp"
#include "matrix_int.hpp"


static inline f64_t rand_f64_t(void) {
    return((rand() % 1000) / 1000.0);
}



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

void matrix::init_rand() {
	for (uint32_t i = 0; i < (this->n * this->m); i++) {
		this->data[i] = rand_f64_t();
	}
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

void matrix::reduce(void) {
	matrix_REF(this->data, this->n, this->m);
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


sq_matrix& sq_matrix::inverse(void) {
	if (!this->_inverse) { this->_inverse = new sq_matrix(this->n); }
	if (!this->RRWS) { this->new_RRWS(); }
	matrix_RREF(this->RRWS->data, this->RRWS->n, this->RRWS->m);
	for (uint32_t i = 0; i < this->RRWS->n; i++) {
		memcpy(&this->_inverse->data[i*m], &this->RRWS->data[i*2*m+m], sizeof(f64_t) * this->m);
	}
	return *this->_inverse;
}


void sq_matrix::LU_decomposition(void) {
	// TODO
}

void sq_matrix::new_RRWS(void) {
	// TODO: CHANGE!!! (use compose.asm)
	RRWS = new matrix(this->n, this->m * 2);
	for (uint32_t i = 0; i < n; i++) {
		memcpy(&RRWS->data[i*m*2], &this->data[i * m], sizeof(f64_t) * m);
		RRWS->data[i*m*2 + m + i] = 1.0;	// identity matrix
	}
}

}