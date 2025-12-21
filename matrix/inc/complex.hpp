//
// Created by marijn on 12/21/25.
//

#ifndef COMPLEX_HPP
#define COMPLEX_HPP
#include "matrix.hpp"



namespace MAT {
	template<class elem_t>
	class complex: public matrix<elem_t, 2, 2> {
	public:
		complex() : matrix<elem_t, 2, 2>() {}
		complex(elem_t a, elem_t b);

		complex& operator=(const complex& rhs);
		complex operator+(const complex& rhs) const;
		complex& operator+=(const complex& rhs);
		complex operator-(const complex& rhs) const;
		complex& operator-=(const complex& rhs);
		complex operator*(const elem_t scalar) const;
		complex& operator*=(const elem_t scalar);
		complex operator*(const complex& rhs) const;
		complex& operator*=(const complex& rhs);

		void print() const;
	};
}




#include <immintrin.h>
namespace MAT {
	template<class elem_t>
	complex<elem_t>::complex(elem_t a, elem_t b) : matrix<elem_t, 2, 2>() {
		this->data[0] =		this->data[3] = a;
		this->data[1] = -b;	this->data[2] = b;
	}

	template<class elem_t>
	complex<elem_t>& complex<elem_t>::operator=(const complex<elem_t>& rhs) {
		memcpy(this->data, rhs.data, 2 * sizeof(elem_t));
		return *this;
	}

	template<class elem_t>
	complex<elem_t> complex<elem_t>::operator+(const complex& rhs) const {
		complex<elem_t> result;
		if constexpr (tf32<elem_t>) {
			v128s_t A = _mm_loadu_ps(this->data);
			const v128s_t B = _mm_loadu_ps(rhs.data);
			A = _mm_add_ps(A, B);	// A += B
			_mm_storeu_ps(result.data, A);
		}
		else if constexpr (tf64<elem_t>) {
			v256d_t A = _mm256_loadu_pd(this->data);
			const v256d_t B = _mm256_loadu_pd(rhs.data);
			A = _mm256_add_pd(A, B);	// A += B
			_mm256_storeu_pd(result.data, A);
		}
		else {
			// TODO
		}
		return result;
	}

	template<class elem_t>
	complex<elem_t>& complex<elem_t>::operator+=(const complex& rhs) {
		if constexpr (tf32<elem_t>) {
			v128s_t A = _mm_loadu_ps(this->data);
			const v128s_t B = _mm_loadu_ps(rhs.data);
			A = _mm_add_ps(A, B);	// A += B
			_mm_storeu_ps(this->data, A);
		}
		else if constexpr (tf64<elem_t>) {
			v256d_t A = _mm256_loadu_pd(this->data);
			const v256d_t B = _mm256_loadu_pd(rhs.data);
			A = _mm256_add_pd(A, B);	// A += B
			_mm256_storeu_pd(this->data, A);
		}
		else {
			// TODO
		}
		return *this;
	}

	template<class elem_t>
	complex<elem_t> complex<elem_t>::operator-(const complex& rhs) const {
		complex<elem_t> result;
		if constexpr (tf32<elem_t>) {
			v128s_t A = _mm_loadu_ps(this->data);
			const v128s_t B = _mm_loadu_ps(rhs.data);
			A = _mm_sub_ps(A, B);	// A += B
			_mm_storeu_ps(result.data, A);
		}
		else if constexpr (tf64<elem_t>) {
			v256d_t A = _mm256_loadu_pd(this->data);
			const v256d_t B = _mm256_loadu_pd(rhs.data);
			A = _mm256_sub_pd(A, B);	// A += B
			_mm256_storeu_pd(result.data, A);
		}
		else {
			// TODO
		}
		return result;
	}

	template<class elem_t>
	complex<elem_t>& complex<elem_t>::operator-=(const complex& rhs) {
		if constexpr (tf32<elem_t>) {
			v128s_t A = _mm_loadu_ps(this->data);
			const v128s_t B = _mm_loadu_ps(rhs.data);
			A = _mm_add_ps(A, B);	// A += B
			_mm_storeu_ps(this->data, A);
		}
		else if constexpr (tf64<elem_t>) {
			v256d_t A = _mm256_loadu_pd(this->data);
			const v256d_t B = _mm256_loadu_pd(rhs.data);
			A = _mm256_add_pd(A, B);	// A += B
			_mm256_storeu_pd(this->data, A);
		}
		else {
			// TODO
		}
		return *this;
	}

	template<class elem_t>
	complex<elem_t> complex<elem_t>::operator*(const elem_t scalar) const {
		complex<elem_t> result;
		if constexpr (tf32<elem_t>) {
			const v128s_t B = _mm_broadcast_ss(&scalar);	// B = [scalar, ..., scalar]
			v128s_t A = _mm_loadu_ps(this->data);
			A = _mm_mul_ps(A, B);	// A *= B
			_mm_storeu_ps(result.data, A);
		}
		else if constexpr (tf64<elem_t>) {
			const v256d_t B = _mm256_broadcast_sd(&scalar);	// B = [scalar, ..., scalar]
			v256d_t A = _mm256_loadu_pd(this->data);
			A = _mm256_mul_pd(A, B);	// A *= B
			_mm256_storeu_pd(result.data, A);
		}
		else {
			// TODO
		}
		return result;
	}

	template<class elem_t>
	complex<elem_t>& complex<elem_t>::operator*=(const elem_t scalar) {
		if constexpr (tf32<elem_t>) {
			const v128s_t B = _mm_broadcast_ss(&scalar);	// B = [scalar, ..., scalar]
			v128s_t A = _mm_loadu_ps(this->data);
			A = _mm_mul_ps(A, B);	// A *= B
			_mm_storeu_ps(this->data, A);
		}
		else if constexpr (tf64<elem_t>) {
			const v256d_t B = _mm256_broadcast_sd(&scalar);	// B = [scalar, ..., scalar]
			v256d_t A = _mm256_loadu_pd(this->data);
			A = _mm256_mul_pd(A, B);	// A *= B
			_mm256_storeu_pd(this->data, A);
		}
		else {
			// TODO
		}
		return *this;
	}

	template<class elem_t>
	complex<elem_t> complex<elem_t>::operator*(const complex& rhs) const {
		complex<elem_t> result;
		if constexpr (tf32<elem_t>) {
			// v128s_t A = _mm_set_ps(this->data[2],	this->data[3]);		// [b, a]
			// v128s_t B = _mm_set_ps(rhs.data[2],		rhs.data[3]);		// [d, c]
			// // A * B = (ac - bd) + i(ad + bc)
			// // [(ac - bd), -(ad + bc)]
			// // [(ad + bc),  (ac - bd)]

			// v128s_t bd_ac = _mm_mul_ps(A, B);		// [bd, ac]
			// B = _mm_shuffle_ps(B, B, 1);	// B = [c, d]
			// v128s_t bc_ad = _mm_mul_ps(A, B);		// [bc, ad]
		}
		else if constexpr (tf64<elem_t>) {
			// TODO
		}
		else {
			// TODO
		}
		return result;
	}

	template<class elem_t>
	complex<elem_t>& complex<elem_t>::operator*=(const complex& rhs) {
		// TODO
		return *this;
	}

	template<class elem_t>
	void complex<elem_t>::print() const {
		printf("%5.2f %c %4.2fi\n",
			this->data[3],
			(this->data[2] > 0) ? '+' : '-',
			(this->data[2] > 0) ? (this->data[2]) : (-this->data[2])
		);
	}
};




#endif //COMPLEX_HPP
