#include "matrix.hpp"

#include <stdexcept>

#include "matrix_int.hpp"


static inline f64_t rand_f64_t(void) {
    return((rand() % 1000) / 1000.0);
}
static inline f32_t rand_f32_t(void) {
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
template<> matrix<f64_t>::matrix(const uint32_t n, const uint32_t m) {
	this->n = n;
	this->m = m;
	this->data = (f64_t*)malloc(sizeof(f64_t) * n * m);
	// TODO: aligned_alloc
}
template<> matrix<f32_t>::matrix(const uint32_t n, const uint32_t m) {
	this->n = n;
	this->m = m;
	this->data = (f32_t*)malloc(sizeof(f32_t) * n * m);
	// TODO: aligned_alloc
}

template<> void matrix<f64_t>::init_rand() {
	for (uint32_t i = 0; i < (this->n * this->m); i++) {
		this->data[i] = rand_f64_t();
	}
}
template<> void matrix<f32_t>::init_rand() {
	for (uint32_t i = 0; i < (this->n * this->m); i++) {
		this->data[i] = rand_f32_t();
	}
}

template<> matrix<f64_t>& matrix<f64_t>::operator=(const matrix<f64_t>& rhs) {
	if (this->m != rhs.m || this->n != rhs.n) { throw std::runtime_error("matrix assignment error"); }
	if (this == &rhs) { return *this; }
	memcpy(this->data, rhs.data, sizeof(f64_t) * this->n * this->m);
	return *this;
}
template<> matrix<f32_t>& matrix<f32_t>::operator=(const matrix<f32_t>& rhs) {
	if (this->m != rhs.m || this->n != rhs.n) { throw std::runtime_error("matrix assignment error"); }
	if (this == &rhs) { return *this; }
	memcpy(this->data, rhs.data, sizeof(f32_t) * this->n * this->m);
	return *this;
}

template<> matrix<f64_t> matrix<f64_t>::operator+(const matrix<f64_t>& rhs) {
	if (this->m != rhs.m || this->n != rhs.n) { throw std::runtime_error("matrix addition error"); }
	matrix<f64_t> result(this->n, this->m);
	matrix_add_64(result.data, this->data, rhs.data, this->n * this->m);
	return result;
}
template<> matrix<f32_t> matrix<f32_t>::operator+(const matrix<f32_t>& rhs) {
	if (this->m != rhs.m || this->n != rhs.n) { throw std::runtime_error("matrix addition error"); }
	matrix<f32_t> result(this->n, this->m);
	matrix_add_32(result.data, this->data, rhs.data, this->n * this->m);
	return result;
}

template<> matrix<f64_t>& matrix<f64_t>::operator+=(const matrix<f64_t>& rhs) {
	if (this->m != rhs.m || this->n != rhs.n) { throw std::runtime_error("matrix addition assignment error"); }
	matrix_add_64(this->data, this->data, rhs.data, this->n * this->m);
	return *this;
}
template<> matrix<f32_t>& matrix<f32_t>::operator+=(const matrix<f32_t>& rhs) {
	if (this->m != rhs.m || this->n != rhs.n) { throw std::runtime_error("matrix addition assignment error"); }
	matrix_add_32(this->data, this->data, rhs.data, this->n * this->m);
	return *this;
}

template<> matrix<f64_t> matrix<f64_t>::operator-(const matrix<f64_t>& rhs) {
	if (this->m != rhs.m || this->n != rhs.n) { throw std::runtime_error("matrix subtraction error"); }
	matrix<f64_t> result(this->n, this->m);
	matrix_sub_64(result.data, this->data, rhs.data, this->n * this->m);
	return result;
}
template<> matrix<f32_t> matrix<f32_t>::operator-(const matrix<f32_t>& rhs) {
	if (this->m != rhs.m || this->n != rhs.n) { throw std::runtime_error("matrix subtraction error"); }
	matrix<f32_t> result(this->n, this->m);
	matrix_sub_32(result.data, this->data, rhs.data, this->n * this->m);
	return result;
}

template<> matrix<f64_t>& matrix<f64_t>::operator-=(const matrix<f64_t>& rhs) {
	if (this->m != rhs.m || this->n != rhs.n) { throw std::runtime_error("matrix subtraction assignment error"); }
	matrix_sub_64(this->data, this->data, rhs.data, this->n * this->m);
	return *this;
}
template<> matrix<f32_t>& matrix<f32_t>::operator-=(const matrix<f32_t>& rhs) {
	if (this->m != rhs.m || this->n != rhs.n) { throw std::runtime_error("matrix subtraction assignment error"); }
	matrix_sub_32(this->data, this->data, rhs.data, this->n * this->m);
	return *this;
}

template<> matrix<f64_t> matrix<f64_t>::operator*(const matrix<f64_t>& rhs) {
	if (this->m != rhs.n) { throw std::runtime_error("matrix multiplication error"); }
	// TODO
	return *this;
}
template<> matrix<f32_t> matrix<f32_t>::operator*(const matrix<f32_t>& rhs) {
	if (this->m != rhs.n) { throw std::runtime_error("matrix multiplication error"); }
	// TODO
	return *this;
}

template<> matrix<f64_t>& matrix<f64_t>::operator*=(const matrix<f64_t>& rhs) {
	if (this->m != rhs.n) { throw std::runtime_error("matrix multiplication assignment error"); }
	// TODO
	return *this;
}
template<> matrix<f32_t>& matrix<f32_t>::operator*=(const matrix<f32_t>& rhs) {
	if (this->m != rhs.n) { throw std::runtime_error("matrix multiplication assignment error"); }
	// TODO
	return *this;
}


template<> matrix<f64_t> matrix<f64_t>::operator*(const f64_t scalar) {
	matrix<f64_t> result(this->n, this->m);
	matrix_scale_64(result.data, this->data, scalar, this->n * this->m);
	return result;
}
template<> matrix<f32_t> matrix<f32_t>::operator*(const f32_t scalar) {
	matrix<f32_t> result(this->n, this->m);
	matrix_scale_32(result.data, this->data, scalar, this->n * this->m);
	return result;
}

template<> matrix<f64_t>& matrix<f64_t>::operator*=(const f64_t scalar) {
	matrix_scale_64(this->data, this->data, scalar, this->n * this->m);
	return *this;
}
template<> matrix<f32_t>& matrix<f32_t>::operator*=(const f32_t scalar) {
	matrix_scale_32(this->data, this->data, scalar, this->n * this->m);
	return *this;

}



template<typename elem_t>
void matrix<elem_t>::print(void) const {
	for (uint32_t i = 0; i < this->n; i++) {
		printf("[");
		for (uint32_t j = 0; j < this->m; j++) {
			printf("%5.2f ", this->data[i * this->m + j]);
		}
		printf("]\n");
	}
	printf("\n");
}
template void matrix<f64_t>::print(void) const;
template void matrix<f32_t>::print(void) const;


template<> void matrix<f64_t>::REF(void) {
	matrix_REF_64(this->data, this->n, this->m);
}
template<> void matrix<f32_t>::REF(void) {
	matrix_REF_32(this->data, this->n, this->m);
}

template<> void matrix<f64_t>::RREF(void) {
	matrix_RREF_64(this->data, this->n, this->m);
}
template<> void matrix<f32_t>::RREF(void) {
	matrix_RREF_32(this->data, this->n, this->m);
}




/*!<
 * square matrix
 * */
//sq_matrix& sq_matrix::L_matrix(void) {
//	if (this->LU[0] == nullptr) { LU_decomposition(); }
//	return *this->LU[0];
//}
//
//sq_matrix& sq_matrix::U_matrix(void) {
//	if (this->LU[1] == nullptr) { LU_decomposition(); }
//	return *this->LU[1];
//}

template<> void sq_matrix<f64_t>::new_RRWS(void) {
	// TODO: CHANGE!!! (use compose.asm)
	RRWS = new matrix<f64_t>(this->n, this->m * 2);
	for (uint32_t i = 0; i < n; i++) {
		memcpy(&RRWS->data[i*m*2], &this->data[i * m], sizeof(f64_t) * m);
		RRWS->data[i*m*2 + m + i] = 1.0;	// identity matrix
	}
}
template<> void sq_matrix<f32_t>::new_RRWS(void) {
	// TODO: CHANGE!!! (use compose.asm)
	RRWS = new matrix<f32_t>(this->n, this->m * 2);
	for (uint32_t i = 0; i < n; i++) {
		memcpy(&RRWS->data[i*m*2], &this->data[i * m], sizeof(f32_t) * m);
		RRWS->data[i*m*2 + m + i] = 1.0;	// identity matrix
	}
}

template<> sq_matrix<f64_t>& sq_matrix<f64_t>::inverse(void) {
	if (!this->_inverse) { this->_inverse = new sq_matrix(this->n); }
	if (!this->RRWS) { this->new_RRWS(); }
	matrix_RREF_64(this->RRWS->data, this->RRWS->n, this->RRWS->m);
	for (uint32_t i = 0; i < this->RRWS->n; i++) {
		memcpy(&this->_inverse->data[i*m], &this->RRWS->data[i*2*m+m], sizeof(f64_t) * this->m);
	}
	return *this->_inverse;
}
template<> sq_matrix<f32_t>& sq_matrix<f32_t>::inverse(void) {
	if (!this->_inverse) { this->_inverse = new sq_matrix(this->n); }
	if (!this->RRWS) { this->new_RRWS(); }
	matrix_RREF_32(this->RRWS->data, this->RRWS->n, this->RRWS->m);
	for (uint32_t i = 0; i < this->RRWS->n; i++) {
		memcpy(&this->_inverse->data[i*m], &this->RRWS->data[i*2*m+m], sizeof(f32_t) * this->m);
	}
	return *this->_inverse;
}



}