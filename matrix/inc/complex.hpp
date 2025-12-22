//
// Created by marijn on 12/21/25.
//

#ifndef COMPLEX_HPP
#define COMPLEX_HPP
#include "matrix.hpp"



namespace MAT {
	template<class elem_t>
	class complex {
	public:
		complex(void) = default;
		complex(elem_t real, elem_t imag);
		//complex(const complex<elem_t>& complex);
		complex(const sq_matrix<elem_t, 2>& mat);

		complex& operator=(const complex& rhs);
		complex operator+(const complex& rhs) const;
		complex& operator+=(const complex& rhs);
		complex operator-(const complex& rhs) const;
		complex& operator-=(const complex& rhs);
		complex operator*(const elem_t scalar) const;
		complex& operator*=(const elem_t scalar);
		complex operator*(const complex& rhs) const;
		complex& operator*=(const complex& rhs);

		sq_matrix<elem_t, 2> mat(void) const;

		void print(void) const;
	private:
		elem_t real;
		elem_t imag;
	};
}




namespace MAT {
	template<class elem_t>
	complex<elem_t>::complex(elem_t real, elem_t imag) {
		this->real = real;
		this->imag = imag;
	}

	template<class elem_t>
	complex<elem_t>::complex(const sq_matrix<elem_t, 2>& mat) {
		this->real = mat.data[3];
		this->imag = mat.data[2];
	}

	template<class elem_t>
	complex<elem_t>& complex<elem_t>::operator=(const complex<elem_t>& rhs) {
		this->real = rhs.real;
		this->imag = rhs.imag;
		return *this;
	}

	template<class elem_t>
	complex<elem_t> complex<elem_t>::operator+(const complex& rhs) const {
		complex<elem_t> result;
		if constexpr (tf32<elem_t>) {
			result.real = this->real + rhs.real;
			result.imag = this->imag + rhs.imag;
		}
		else if constexpr (tf64<elem_t>) {
			v128d_t A = _mm_loadu_pd(&this->real);
			const v128d_t B = _mm_loadu_pd(&rhs.real);
			A = _mm_add_pd(A, B);	// A += B
			_mm_storeu_pd(&result.real, A);
		}
		else {
			// TODO
		}
		return result;
	}

	template<class elem_t>
	complex<elem_t>& complex<elem_t>::operator+=(const complex& rhs) {
		if constexpr (tf32<elem_t>) {
			this->real += rhs.real;
			this->imag += rhs.imag;
		}
		else if constexpr (tf64<elem_t>) {
			v128d_t A = _mm_loadu_pd(&this->real);
			const v128d_t B = _mm_loadu_pd(&rhs.real);
			A = _mm_add_pd(A, B);	// A += B
			_mm_storeu_pd(&this->real, A);
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
			result.real = this->real - rhs.real;
			result.imag = this->imag - rhs.imag;
		}
		else if constexpr (tf64<elem_t>) {
			v128d_t A = _mm_loadu_pd(&this->real);
			const v128d_t B = _mm_loadu_pd(&rhs.real);
			A = _mm_sub_pd(A, B);	// A -= B
			_mm_storeu_pd(&result.real, A);
		}
		else {
			// TODO
		}
		return result;
	}

	template<class elem_t>
	complex<elem_t>& complex<elem_t>::operator-=(const complex& rhs) {
		if constexpr (tf32<elem_t>) {
			this->real -= rhs.real;
			this->imag -= rhs.imag;
		}
		else if constexpr (tf64<elem_t>) {
			v128d_t A = _mm_loadu_pd(&this->real);
			const v128d_t B = _mm_loadu_pd(&rhs.real);
			A = _mm_sub_pd(A, B);	// A -= B
			_mm_storeu_pd(&this->real, A);
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
			result.real = this->real * scalar;
			result.imag = this->imag * scalar;
		}
		else if constexpr (tf64<elem_t>) {
			const v128d_t B = _mm_set_pd(scalar, scalar);	// B = [scalar, scalar]
			v128d_t A = _mm_loadu_pd(&this->real);
			A = _mm_mul_pd(A, B);	// A *= B
			_mm_storeu_pd(&result.real, A);
		}
		else {
			// TODO
		}
		return result;
	}

	template<class elem_t>
	complex<elem_t>& complex<elem_t>::operator*=(const elem_t scalar) {
		if constexpr (tf32<elem_t>) {
			this->real *= scalar;
			this->imag *= scalar;
		}
		else if constexpr (tf64<elem_t>) {
			const v128d_t B = _mm_set_pd(scalar, scalar);	// B = [scalar, scalar]
			v128d_t A = _mm_loadu_pd(&this->real);
			A = _mm_mul_pd(A, B);	// A *= B
			_mm_storeu_pd(&this->real, A);
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
			v128s_t A =		_mm_set_ps(0.0f, 0.0f, this->imag,	this->real);		// [0, 0, imag (b), real (a)]
			v128s_t B =		_mm_set_ps(0.0f, 0.0f, rhs.imag,	rhs.real);		// [0, 0, imag (d), real (c)]
			v128s_t M1 =	_mm_mul_ps(A, B);		// [bd, ac]
			v128s_t M2 =	_mm_shuffle_ps(B, B, _MM_SHUFFLE(2,3,0,1));
			M2 =			_mm_mul_ps(A, M2);	// [bc, ad]
			M1 = _mm_sub_ss(	// ac - bd
				M1,
				_mm_shuffle_ps(M1, M1, _MM_SHUFFLE(2,3,0,1))
			);
			M2 = _mm_add_ss(	// ad + bc
				M2,
				_mm_shuffle_ps(M2, M2, _MM_SHUFFLE(2,3,0,1))
			);

			result.real = _mm_cvtss_f32(M1);
			result.imag = _mm_cvtss_f32(M2);
		}
		else if constexpr (tf64<elem_t>) {
			v128d_t A =		_mm_set_pd(this->imag,	this->real);	// [imag (b), real (a)]
			v128d_t B =		_mm_set_pd(rhs.imag,		rhs.real);	// [imag (d), real (c)]
			v128d_t M1 =	_mm_mul_pd(A, B);		// [bd, ac]
			v128d_t M2 =	_mm_shuffle_pd(B, B, 1);
			M2 =			_mm_mul_pd(A, M2);	// [bc, ad]
			M1 = _mm_sub_sd(	// ac - bd
				M1,
				_mm_shuffle_pd(M1, M1, 1)
			);
			M2 = _mm_add_sd(	// ad + bc
				M2,
				_mm_shuffle_pd(M2, M2, 1)
			);

			result.real = _mm_cvtsd_f64(M1);
			result.imag = _mm_cvtsd_f64(M2);
		}
		else {
			// TODO
		}
		return result;
	}

	template<class elem_t>
	complex<elem_t>& complex<elem_t>::operator*=(const complex& rhs) {
				if constexpr (tf32<elem_t>) {
			v128s_t A =		_mm_set_ps(0.0f, 0.0f, this->imag,	this->real);		// [0, 0, imag (b), real (a)]
			v128s_t B =		_mm_set_ps(0.0f, 0.0f, rhs.imag,	rhs.real);		// [0, 0, imag (d), real (c)]
			v128s_t M1 =	_mm_mul_ps(A, B);		// [bd, ac]
			v128s_t M2 =	_mm_shuffle_ps(B, B, _MM_SHUFFLE(2,3,0,1));
			M2 =			_mm_mul_ps(A, M2);	// [bc, ad]
			M1 = _mm_sub_ss(	// ac - bd
				M1,
				_mm_shuffle_ps(M1, M1, _MM_SHUFFLE(2,3,0,1))
			);
			M2 = _mm_add_ss(	// ad + bc
				M2,
				_mm_shuffle_ps(M2, M2, _MM_SHUFFLE(2,3,0,1))
			);

			this->real = _mm_cvtss_f32(M1);
			this->imag = _mm_cvtss_f32(M2);
		}
		else if constexpr (tf64<elem_t>) {
			v128d_t A =		_mm_set_pd(this->imag,	this->real);	// [imag (b), real (a)]
			v128d_t B =		_mm_set_pd(rhs.imag,		rhs.real);	// [imag (d), real (c)]
			v128d_t M1 =	_mm_mul_pd(A, B);		// [bd, ac]
			v128d_t M2 =	_mm_shuffle_pd(B, B, 1);
			M2 =			_mm_mul_pd(A, M2);	// [bc, ad]
			M1 = _mm_sub_sd(	// ac - bd
				M1,
				_mm_shuffle_pd(M1, M1, 1)
			);
			M2 = _mm_add_sd(	// ad + bc
				M2,
				_mm_shuffle_pd(M2, M2, 1)
			);

			this->real = _mm_cvtsd_f64(M1);
			this->imag = _mm_cvtsd_f64(M2);
		}
		else {
			// TODO
		}
		return *this;
	}

	template<class elem_t>
	sq_matrix<elem_t, 2> complex<elem_t>::mat(void) const {
		sq_matrix<elem_t, 2> result;
		result.data[0] = result.data[3] = this->real;
		result.data[1] = -this->imag;
		result.data[2] = this->imag;
		return result;
	}

	template<class elem_t>
	void complex<elem_t>::print(void) const {
		printf("%5.2f %c %4.2fi\n",
			this->real,
			(this->imag > 0) ? '+' : '-',
			(this->imag > 0) ? (this->imag) : (-this->imag)
		);
	}
};




#endif //COMPLEX_HPP
