#ifndef MATRIX_LIBRARY_H
#define MATRIX_LIBRARY_H
#include <stdint.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "types.h"


namespace MAT {
class sq_matrix;


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
class matrix {
friend class sq_matrix;
public:
	// TODO: copy constructor: matrix(const matrix&) = default;
	matrix(uint32_t n, uint32_t m);
	~matrix(void);

	void init_rand(void);

	void REF(void);
	void REF_C(void);
	void RREF(void);

	f64_t& operator()(uint32_t i, uint32_t j);
	void print(void) const;

private:
	uint32_t n, m;
	f64_t* data;
};


/*!<
 * square matrix
 * */
class sq_matrix : public matrix {
public:
	sq_matrix(const uint32_t dim) : matrix(dim, dim) {};

	sq_matrix& L_matrix(void);
	sq_matrix& U_matrix(void);
	matrix& get_RRWS(void) { if (!this->RRWS) { this->new_RRWS(); } return *this->RRWS; }  // TODO: this is purely for debug purposes

	sq_matrix& inverse(void);
	//void row_reduce(void);

private:
	void LU_decomposition(void);
	void new_RRWS(void);

	det_t det{};
	sq_matrix* LU[2] = {nullptr, nullptr};
	sq_matrix* _inverse = nullptr;
	matrix* RRWS = nullptr;  // row reduction workspace
};

}



#endif //MATRIX_LIBRARY_H