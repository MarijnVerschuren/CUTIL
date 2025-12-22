//
// Created by marijn on 12/17/25.
//

#ifndef SQ_MATRIX_HPP
#define SQ_MATRIX_HPP
#include "matrix.hpp"


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
	 * square matrix
	 * */
	template<class elem_t, uint32_t n>
	class sq_matrix : public matrix<elem_t, n, n> {
	public:
		sq_matrix(void);
		~sq_matrix(void) = default;

		//sq_matrix& inverse(void);

	private:
		//void LU_decomposition(void);
		//void new_RRWS(void);

		det_t det{};
		//sq_matrix* LU[2] = {nullptr, nullptr};
		//sq_matrix* _inverse = nullptr;
		//matrix<elem_t>* RRWS = nullptr;  // row reduction workspace
	};
}


namespace MAT {
	/*!<
	 * Property classes
	 * */
	inline f64_t& det_t::operator =(const f64_t val) {
		return this->_det = val;
	}




	/*!<
	 * square matrix
	 * */
	template<class elem_t, uint32_t n>
	sq_matrix<elem_t, n>::sq_matrix(void) : matrix<elem_t, n, n>() {}

	//sq_matrix& sq_matrix::L_matrix(void) {
	//	if (this->LU[0] == nullptr) { LU_decomposition(); }
	//	return *this->LU[0];
	//}
	//
	//sq_matrix& sq_matrix::U_matrix(void) {
	//	if (this->LU[1] == nullptr) { LU_decomposition(); }
	//	return *this->LU[1];
	//}

	// template<class elem_t>  void sq_matrix<elem_t>::new_RRWS(void) {
	// 	// TODO: CHANGE!!! (use compose.asm)
	// 	RRWS = new matrix<elem_t>(this->n, this->m * 2);
	// 	for (uint32_t i = 0; i < this->n; i++) {
	// 		memcpy(&RRWS->data[i*this->m*2], &this->data[i * this->m], sizeof(elem_t) * this->m);
	// 		RRWS->data[i*this->m*2 + this->m + i] = 1.0;	// identity matrix
	// 	}
	// }
	//
	// template<> inline sq_matrix<f64_t>& sq_matrix<f64_t>::inverse(void) {
	// 	if (!this->_inverse) { this->_inverse = new sq_matrix(this->n); }
	// 	if (!this->RRWS) { this->new_RRWS(); }
	// 	matrix_RREF_64(this->RRWS->data, this->RRWS->n, this->RRWS->m);
	// 	for (uint32_t i = 0; i < this->RRWS->n; i++) {
	// 		memcpy(&this->_inverse->data[i*m], &this->RRWS->data[i*2*m+m], sizeof(f64_t) * this->m);
	// 	}
	// 	return *this->_inverse;
	// }
	// template<> inline sq_matrix<f32_t>& sq_matrix<f32_t>::inverse(void) {
	// 	if (!this->_inverse) { this->_inverse = new sq_matrix(this->n); }
	// 	if (!this->RRWS) { this->new_RRWS(); }
	// 	matrix_RREF_32(this->RRWS->data, this->RRWS->n, this->RRWS->m);
	// 	for (uint32_t i = 0; i < this->RRWS->n; i++) {
	// 		memcpy(&this->_inverse->data[i*m], &this->RRWS->data[i*2*m+m], sizeof(f32_t) * this->m);
	// 	}
	// 	return *this->_inverse;
	// }
}


#endif //SQ_MATRIX_HPP
