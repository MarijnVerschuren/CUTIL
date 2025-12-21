#ifndef MATRIX_LIBRARY_H
#define MATRIX_LIBRARY_H
#include <stdint.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "types.h"
#include "matrix_int.hpp"



template<typename T>
concept tf32 = std::is_same_v<T, f32_t>;
template<typename T>
concept tf64 = std::is_same_v<T, f64_t>;

#define RAND() ((rand() % 1000) / 1000.0)




namespace MAT {
	/*!<
	 * matrix
	 * */
	template<class elem_t, u32_t m, u32_t n>
	class matrix {
	template <class, u32_t> friend class sq_matrix;
	template <class, u32_t, u32_t> friend class matrix;
	template <class> friend class complex;
	public:
		// TODO: copy constructor: matrix(const matrix&) = default;
		matrix();
		~matrix(void) { free(this->data); }

		void init_rand(void);

		void REF(void);
		void RREF(void);

		elem_t& operator()(uint32_t i, uint32_t j) { return this->data[i * m + j]; }
		matrix& operator=(const matrix& rhs);
		matrix operator+(const matrix& rhs) const;
		matrix& operator+=(const matrix& rhs);
		matrix operator-(const matrix& rhs) const;
		matrix& operator-=(const matrix& rhs);
		matrix operator*(const elem_t scalar) const;	// scale
		matrix& operator*=(const elem_t scalar);		// scale

		template<u32_t p>
		matrix<elem_t, m, p> operator*(const matrix<elem_t, n, p>& rhs) const;

		void print(void) const;

	private:
		//uint8_t* null;	// bit array
		elem_t* data;
	};
}



namespace MAT {
	/*!<
	 * matrix
	 * */
	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n>::matrix() {
		this->data = (elem_t*)malloc(sizeof(elem_t) * n * m);
	}


	template<class elem_t, u32_t m, u32_t n>
	void matrix<elem_t, m, n>::init_rand() {
		for (uint32_t i = 0; i < (n * m); i++) {
			this->data[i] = RAND();
		}
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n>& matrix<elem_t, m, n>::operator=(const matrix<elem_t, m, n>& rhs) {
		if (this == &rhs) { return *this; }
		memcpy(this->data, rhs.data, sizeof(elem_t) * n * m);
		return *this;
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n> matrix<elem_t, m, n>::operator+(const matrix<elem_t, m, n>& rhs) const {
		matrix<elem_t, m, n> result;
		if constexpr (tf32<elem_t>)			{ matrix_add_32(result.data, this->data, rhs.data, n * m); }
		else if constexpr (tf64<elem_t>)	{ matrix_add_64(result.data, this->data, rhs.data, n * m); }
		else {
			// TODO: C algo!
		}
		return result;
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n>& matrix<elem_t, m, n>::operator+=(const matrix<elem_t, m, n>& rhs) {
		if constexpr (tf32<elem_t>)			{ matrix_add_32(this->data, this->data, rhs.data, n * m); }
		else if constexpr (tf64<elem_t>)	{ matrix_add_64(this->data, this->data, rhs.data, n * m); }
		else {
			// TODO: C algo!
		}
		return *this;
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n> matrix<elem_t, m, n>::operator-(const matrix<elem_t, m, n>& rhs) const {
		matrix<elem_t, m, n> result;
		if constexpr (tf32<elem_t>)			{ matrix_sub_32(result.data, this->data, rhs.data, n * m); }
		else if constexpr (tf64<elem_t>)	{ matrix_sub_64(result.data, this->data, rhs.data, n * m); }
		else {
			// TODO: C algo!
		}
		return result;
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n>& matrix<elem_t, m, n>::operator-=(const matrix<elem_t, m, n>& rhs) {
		if constexpr (tf32<elem_t>)			{ matrix_sub_32(this->data, this->data, rhs.data, n * m); }
		else if constexpr (tf64<elem_t>)	{ matrix_sub_64(this->data, this->data, rhs.data, n * m); }
		else {
			// TODO: C algo!
		}
		return *this;
	}

	template<class elem_t, u32_t m, u32_t n>
	template<u32_t p>
	matrix<elem_t, m, p> matrix<elem_t, m, n>::operator*(const matrix<elem_t, n, p>& rhs) const {
	    matrix<elem_t, m, p> result;

		// TODO: improved LASER method
		// TODO: alman williams matrix multiplication!!!
		if constexpr (tf32<elem_t>)			{ matrix_mul_32(result.data, this->data, rhs.data, m, n, p); }
		else if constexpr (tf64<elem_t>)	{ matrix_mul_64(result.data, this->data, rhs.data, m, n, p); }
		else {
			// TODO: general algo!
		}

	    return result;
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n> matrix<elem_t, m, n>::operator*(const elem_t scalar) const {
		matrix<elem_t, m, n> result;
		if constexpr (tf32<elem_t>)			{ matrix_scale_32(result.data, this->data, scalar, n * m); }
		else if constexpr (tf64<elem_t>)	{ matrix_scale_64(result.data, this->data, scalar, n * m); }
		else {
			// TODO: C algo!
		}
		return result;
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n>& matrix<elem_t, m, n>::operator*=(const elem_t scalar) {
		if constexpr (tf32<elem_t>)			{ matrix_scale_32(this->data, this->data, scalar, n * m); }
		else if constexpr (tf64<elem_t>)	{ matrix_scale_64(this->data, this->data, scalar, n * m); }
		else {
			// TODO: C algo!
		}
		return *this;
	}


	template<class elem_t, u32_t m, u32_t n>
	void matrix<elem_t, m, n>::print(void) const {
		for (uint32_t i = 0; i < m; i++) {
			printf("[");
			for (uint32_t j = 0; j < n; j++) {
				printf("%5.2f ", this->data[i * n + j]);
			}
			printf("]\n");
		}
		printf("\n");
	}


	template<class elem_t, u32_t m, u32_t n>
	void matrix<elem_t, m, n>::REF(void) {
		if constexpr (tf32<elem_t>)			{ matrix_REF_32(this->data, m, n); }
		else if constexpr (tf64<elem_t>)	{ matrix_REF_64(this->data, m, n); }
		else {
			// TODO: C algo!
		}
	}

	template<class elem_t, u32_t m, u32_t n>
	void matrix<elem_t, m, n>::RREF(void) {
		if constexpr (tf32<elem_t>)			{ matrix_RREF_32(this->data, m, n); }
		else if constexpr (tf64<elem_t>)	{ matrix_RREF_64(this->data, m, n); }
		else {
			// TODO: C algo!
		}
	}
}


#endif //MATRIX_LIBRARY_H