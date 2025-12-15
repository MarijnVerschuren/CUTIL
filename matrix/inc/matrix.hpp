#ifndef MATRIX_LIBRARY_H
#define MATRIX_LIBRARY_H
#include <stdint.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "types.h"


namespace MAT {


/*!<
 * Property classes
 * */
class det_t {
public:
	det_t() = default;
	operator f64_t() const { return this->_det; }
	// TODO make framework with which redundant calculations can be avoided (feed values)
	f64_t& operator =(f64_t val);
private:
	f64_t _det;
};


/*!<
 * matrix
 * */
template<typename elem_t>
class matrix {
template <typename> friend class sq_matrix;
public:
	// TODO: copy constructor: matrix(const matrix&) = default;
	matrix(uint32_t n, uint32_t m);
	~matrix(void) { free(this->data); }

	void init_rand(void);

	void REF(void);
	void RREF(void);

	elem_t& operator()(uint32_t i, uint32_t j) { return this->data[i * m + j]; }
	void print(void) const;

private:
	uint32_t n, m;
	elem_t* data;
};
template class matrix<f64_t>;
template class matrix<f32_t>;

/*!<
 * square matrix
 * */
template<typename elem_t>
class sq_matrix : public matrix<elem_t> {
public:
	sq_matrix(const uint32_t dim) : matrix<elem_t>(dim, dim) {};

	//sq_matrix& L_matrix(void);
	//sq_matrix& U_matrix(void);
	matrix<elem_t>& get_RRWS(void) { if (!this->RRWS) { this->new_RRWS(); } return *this->RRWS; }  // TODO: this is purely for debug purposes

	sq_matrix& inverse(void);
	//void row_reduce(void);

private:
	void LU_decomposition(void);
	void new_RRWS(void);

	det_t det{};
	//sq_matrix* LU[2] = {nullptr, nullptr};
	sq_matrix* _inverse = nullptr;
	matrix<elem_t>* RRWS = nullptr;  // row reduction workspace
};
template class sq_matrix<f64_t>;
template class sq_matrix<f32_t>;

}



#endif //MATRIX_LIBRARY_H