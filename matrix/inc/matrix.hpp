#ifndef MATRIX_LIBRARY_H
#define MATRIX_LIBRARY_H
#include <stdint.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "types.h"



template<typename T>
concept tf32 = std::is_same_v<T, f32_t>;
template<typename T>
concept tf64 = std::is_same_v<T, f64_t>;

template<std::size_t m, std::size_t n>
concept sq = (m == n);

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
		static constexpr u8_t big = (m * n * sizeof(elem_t)) > 4096;
		using data_t = std::conditional_t<big, elem_t*, elem_t[m*n]>;

		// TODO: copy constructor: matrix(const matrix&) = default;
		matrix(void);
		~matrix(void);

		void init_rand(void);
		void init_identity(void) requires sq<m, n>;

		void REF(void);
		void RREF(void);

		elem_t& operator()(u32_t i, u32_t j) { return this->data[i * m + j]; }
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
		data_t data;
	};
}



// matrix multiplication block size parameters
#define BM 64
#define BN 64
#define BK 64
namespace MAT {
	/*!<
	 * matrix
	 * */
	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n>::matrix() {
		if constexpr (big) {
			this->data = (elem_t*)malloc(sizeof(elem_t) * n * m);
		}
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n>::~matrix() { if constexpr (big) { free(this->data); } }

	template<class elem_t, u32_t m, u32_t n>
	void matrix<elem_t, m, n>::init_rand() {
		for (u32_t i = 0; i < (n * m); i++) {
			this->data[i] = RAND();
		}
	}

	template<class elem_t, u32_t m, u32_t n>
	void matrix<elem_t, m, n>::init_identity(void) requires sq<m, n> {
		// TODO: (now sq class is not needed!)
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
		u32_t i, j;
		for (i = 0; i < m; i++) {
			j = 0;
			if constexpr (tf32<elem_t> && n > 7) {
				for (; j < n; j += 8) {
					v256s_t A =	_mm256_loadu_ps(&this->data[i * m + j]);
					v256s_t B =	_mm256_loadu_ps(&rhs.data[i * m + j]);
					A =			_mm256_add_ps(A, B);
					_mm256_storeu_ps(&result.data[i * m + j], A);
				}
			} else if constexpr (tf64<elem_t> && n > 3) {
				for (; j < n; j += 4) {
					v256d_t A =	_mm256_loadu_pd(&this->data[i * m + j]);
					v256d_t B =	_mm256_loadu_pd(&rhs.data[i * m + j]);
					A =			_mm256_add_pd(A, B);
					_mm256_storeu_pd(&result.data[i * m + j], A);
				}
			}
			for (; j < n; j++) {
				result.data[i * m + j] = this->data[i * m + j] + rhs.data[i * m + j];
			}
		}
		return result;
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n>& matrix<elem_t, m, n>::operator+=(const matrix<elem_t, m, n>& rhs) {
		u32_t i, j;
		for (i = 0; i < m; i++) {
			j = 0;
			if constexpr (tf32<elem_t> && n > 7) {
				for (; j < n; j += 8) {
					v256s_t A =	_mm256_loadu_ps(&this->data[i * m + j]);
					v256s_t B =	_mm256_loadu_ps(&rhs.data[i * m + j]);
					A =			_mm256_add_ps(A, B);
					_mm256_storeu_ps(&this->data[i * m + j], A);
				}
			} else if constexpr (tf64<elem_t> && n > 3) {
				for (; j < n; j += 4) {
					v256d_t A =	_mm256_loadu_pd(&this->data[i * m + j]);
					v256d_t B =	_mm256_loadu_pd(&rhs.data[i * m + j]);
					A =			_mm256_add_pd(A, B);
					_mm256_storeu_pd(&this->data[i * m + j], A);
				}
			}
			for (; j < n; j++) {
				this->data[i * m + j] += rhs.data[i * m + j];
			}
		}
		return *this;
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n> matrix<elem_t, m, n>::operator-(const matrix<elem_t, m, n>& rhs) const {
		matrix<elem_t, m, n> result;
		u32_t i, j;
		for (i = 0; i < m; i++) {
			j = 0;
			if constexpr (tf32<elem_t> && n > 7) {
				for (; j < n; j += 8) {
					v256s_t A =	_mm256_loadu_ps(&this->data[i * m + j]);
					v256s_t B =	_mm256_loadu_ps(&rhs.data[i * m + j]);
					A =			_mm256_sub_ps(A, B);
					_mm256_storeu_ps(&result.data[i * m + j], A);
				}
			} else if constexpr (tf64<elem_t> && n > 3) {
				for (; j < n; j += 4) {
					v256d_t A =	_mm256_loadu_pd(&this->data[i * m + j]);
					v256d_t B =	_mm256_loadu_pd(&rhs.data[i * m + j]);
					A =			_mm256_sub_pd(A, B);
					_mm256_storeu_pd(&result.data[i * m + j], A);
				}
			}
			for (; j < n; j++) {
				result.data[i * m + j] = this->data[i * m + j] - rhs.data[i * m + j];
			}
		}
		return result;
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n>& matrix<elem_t, m, n>::operator-=(const matrix<elem_t, m, n>& rhs) {
		u32_t i, j;
		for (i = 0; i < m; i++) {
			j = 0;
			if constexpr (tf32<elem_t> && n > 7) {
				for (; j < n; j += 8) {
					v256s_t A =	_mm256_loadu_ps(&this->data[i * m + j]);
					v256s_t B =	_mm256_loadu_ps(&rhs.data[i * m + j]);
					A =			_mm256_add_ps(A, B);
					_mm256_storeu_ps(&this->data[i * m + j], A);
				}
			} else if constexpr (tf64<elem_t> && n > 3) {
				for (; j < n; j += 4) {
					v256d_t A =	_mm256_loadu_pd(&this->data[i * m + j]);
					v256d_t B =	_mm256_loadu_pd(&rhs.data[i * m + j]);
					A =			_mm256_add_pd(A, B);
					_mm256_storeu_pd(&this->data[i * m + j], A);
				}
			}
			for (; j < n; j++) {
				this->data[i * m + j] += rhs.data[i * m + j];
			}
		}
		return *this;
	}

	template<class elem_t, u32_t m, u32_t n>
	template<u32_t p>
	matrix<elem_t, m, p> matrix<elem_t, m, n>::operator*(const matrix<elem_t, n, p>& rhs) const {
	    matrix<elem_t, m, p> result;
		u32_t i, j, k, ii, jj, kk, i_max, j_max, k_max;
		elem_t a, *c_ptr; const elem_t* b_ptr;
		{	// delete zero after zeroing loop
	    	i = 0;
	    	if constexpr (tf32<elem_t>) {
				v256s_t zero = _mm256_setzero_ps();
	    		for (; i + 3 < n*p; i += 4) { _mm256_storeu_ps(&result.data[i], zero); }
	    	} else if constexpr (tf64<elem_t>) {
				v256d_t zero = _mm256_setzero_pd();
	    		for (; i + 3 < n*p; i += 4) { _mm256_storeu_pd(&result.data[i], zero); }
	    	}
	    	for (; i < n*p; ++i) { result.data[i] = 0.0; }
		}

		// block loop
	    for (ii = 0; ii < m; ii += BM) {
		    for (kk = 0; kk < n; kk += BK) {
		    	for (jj = 0; jj < p; jj += BN) {
	                i_max = (ii + BM < m ? ii + BM : m);
	                j_max = (jj + BN < p ? jj + BN : p);
	                k_max = (kk + BK < n ? kk + BK : n);
		    		// TODO: multithread for large matrecies
		    		// task
	                for (i = ii; i < i_max; ++i) {
		                for (k = kk; k < k_max; ++k) {
							a = this->data[i * n + k];
		                	b_ptr = rhs.data + k * p + jj;
							c_ptr = result.data + i * p + jj;

		                	j = 0;
							if constexpr (tf64<elem_t>) {
								v256d_t a_vec, b_vec, c_vec;
								a_vec = _mm256_set1_pd(a); // [a, a, a, a]
								for (; j + 3 < j_max - jj; j += 4) {
									c_vec = _mm256_loadu_pd(c_ptr + j);
									b_vec = _mm256_loadu_pd(b_ptr + j);
									c_vec = _mm256_fmadd_pd(a_vec, b_vec, c_vec); // c += a*b
									_mm256_storeu_pd(c_ptr + j, c_vec);
								}
							} else if constexpr (tf32<elem_t>) {
								v256s_t a_vec, b_vec, c_vec;
								a_vec = _mm256_set1_ps(a); // [a, a, a, a, a, a, a, a]
								for (; j + 7 < j_max - jj; j += 8) {
									c_vec = _mm256_loadu_ps(c_ptr + j);
									b_vec = _mm256_loadu_ps(b_ptr + j);
									c_vec = _mm256_fmadd_ps(a_vec, b_vec, c_vec); // c += a*b
									_mm256_storeu_ps(c_ptr + j, c_vec);
								}
							}
	                        for (; j < j_max - jj; ++j) {
	                            c_ptr[j] += a * b_ptr[j];
	                        }
		                }
	                }
		    		// ~task
		    	}
		    }
	    }

	    return result;
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n> matrix<elem_t, m, n>::operator*(const elem_t scalar) const {
		matrix<elem_t, m, n> result;
		u32_t i, j;
		if constexpr (tf32<elem_t> && n > 7) {
			v256s_t B = _mm256_set1_ps(scalar);
			for (i = 0; i < m; i++) {
				j = 0;
				for (; j < n; j += 8) {
					v256s_t A =	_mm256_loadu_ps(&this->data[i * m + j]);
					A =			_mm256_mul_ps(A, B);
					_mm256_storeu_ps(&result.data[i * m + j], A);
				}
				for (; j < n; j++) {
					result.data[i * m + j] = this->data[i * m + j] * scalar;
				}
			}
		}
		else if constexpr (tf64<elem_t> && n > 3) {
			v256d_t B = _mm256_set1_pd(scalar);
			for (i = 0; i < m; i++) {
				j = 0;
				for (; j < n; j += 4) {
					v256d_t A =	_mm256_loadu_pd(&this->data[i * m + j]);
					A =			_mm256_mul_pd(A, B);
					_mm256_storeu_pd(&result.data[i * m + j], A);
				}
				for (; j < n; j++) {
					result.data[i * m + j] = this->data[i * m + j] * scalar;
				}
			}
		}
		else {
			for (i = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					result.data[i * m + j] = this->data[i * m + j] * scalar;
				}
			}
		}
		return result;
	}

	template<class elem_t, u32_t m, u32_t n>
	matrix<elem_t, m, n>& matrix<elem_t, m, n>::operator*=(const elem_t scalar) {
		u32_t i, j;
		if constexpr (tf32<elem_t> && n > 7) {
			v256s_t B = _mm256_set1_ps(scalar);
			for (i = 0; i < m; i++) {
				j = 0;
				for (; j < n; j += 8) {
					v256s_t A =	_mm256_loadu_ps(&this->data[i * m + j]);
					A =			_mm256_mul_ps(A, B);
					_mm256_storeu_ps(&this->data[i * m + j], A);
				}
				for (; j < n; j++) {
					this->data[i * m + j] *= scalar;
				}
			}
		}
		else if constexpr (tf64<elem_t> && n > 3) {
			v256d_t B = _mm256_set1_pd(scalar);
			for (i = 0; i < m; i++) {
				j = 0;
				for (; j < n; j += 4) {
					v256d_t A =	_mm256_loadu_pd(&this->data[i * m + j]);
					A =			_mm256_mul_pd(A, B);
					_mm256_storeu_pd(&this->data[i * m + j], A);
				}
				for (; j < n; j++) {
					this->data[i * m + j] *= scalar;
				}
			}
		}
		else {
			for (i = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					this->data[i * m + j] *= scalar;
				}
			}
		}
		return *this;
	}

	template<class elem_t, u32_t m, u32_t n>
	void matrix<elem_t, m, n>::print(void) const {
		for (u32_t i = 0; i < m; i++) {
			printf("[");
			for (u32_t j = 0; j < n; j++) {
				printf("%5.2f ", this->data[i * n + j]);
			}
			printf("]\n");
		}
		printf("\n");
	}

	template<class elem_t, u32_t m, u32_t n>
	void matrix<elem_t, m, n>::REF(void) {
		u32_t i, j, k;
		elem_t pivot, tmp;
		for (i = 0; i < m; i++) {
			for (j = i; j < m; j++) {
				pivot = this->data[j * n + i];
				if (pivot != 0.0f) { break; }
			}
			if (pivot == 0.0f) { continue; }
			if (i != j) { // swap rows
				k = 0;
				if constexpr (tf32<elem_t> && n > 7) {
					v256s_t A, B;
					for (k = j; k < n; k += 8) {
						A =	_mm256_loadu_ps(&this->data[i * n + k]);
						B =	_mm256_loadu_ps(&this->data[j * n + k]);
						_mm256_storeu_ps(&this->data[i * n + k], B);
						_mm256_storeu_ps(&this->data[j * n + k], A);
					}
				} else if constexpr (tf64<elem_t> && n > 3) {
					v256d_t A, B;
					for (k = j; k < n; k += 4) {
						A =	_mm256_loadu_pd(&this->data[i * n + k]);
						B =	_mm256_loadu_pd(&this->data[j * n + k]);
						_mm256_storeu_pd(&this->data[i * n + k], B);
						_mm256_storeu_pd(&this->data[j * n + k], A);
					}
				}
				for (; k < n; k++) {
					tmp =					this->data[i * n + k];
					this->data[i * n + k] = this->data[j * n + k];
					this->data[j * n + k] = tmp;
				}
			}
			for (j = i + 1; j < n; j++) {
				k = 0;
				tmp = -this->data[j * n + i]/pivot;
				if constexpr (tf32<elem_t> && n > 7) {
					v256s_t A, B, C = _mm256_set1_ps(tmp);
					for (; k < n; k += 8) {
						A =	_mm256_loadu_ps(&this->data[j * n + k]);
						B =	_mm256_loadu_ps(&this->data[i * n + k]);
						A = _mm256_fmadd_ps(B, C, A);	// B * C + A
						_mm256_storeu_ps(&this->data[j * n + k], A);
					}
				} else if constexpr (tf64<elem_t> && n > 3) {
					v256d_t A, B, C = _mm256_set1_pd(tmp);
					for (; k < n; k += 4) {
						A =	_mm256_loadu_pd(&this->data[j * n + k]);
						B =	_mm256_loadu_pd(&this->data[i * n + k]);
						A = _mm256_fmadd_pd(B, C, A);	// B * C + A
						_mm256_storeu_pd(&this->data[j * n + k], A);
					}
				}
				for (; k < n; k++) {
					this->data[j * m + k] +=  this->data[i * m + k] * tmp;
				}
			}
		}
	}

	template<class elem_t, u32_t m, u32_t n>
	void matrix<elem_t, m, n>::RREF(void) {
		u32_t i, j, k;
		elem_t pivot, tmp;
		for (i = 0; i < m; i++) {
			for (j = i; j < m; j++) {
				pivot = data[j * n + i];
				if (pivot != 0.0f) { break; }
			}
			if (pivot == 0.0f) { continue; }
			if (i != j) { // swap rows
				k = 0;
				if constexpr (tf32<elem_t> && n > 7) {
					v256s_t A, B;
					for (k = j; k < n; k += 8) {
						A =	_mm256_loadu_ps(&this->data[i * n + k]);
						B =	_mm256_loadu_ps(&this->data[j * n + k]);
						_mm256_storeu_ps(&this->data[i * n + k], B);
						_mm256_storeu_ps(&this->data[j * n + k], A);
					}
				} else if constexpr (tf64<elem_t> && n > 3) {
					v256d_t A, B;
					for (k = j; k < n; k += 4) {
						A =	_mm256_loadu_pd(&this->data[i * n + k]);
						B =	_mm256_loadu_pd(&this->data[j * n + k]);
						_mm256_storeu_pd(&this->data[i * n + k], B);
						_mm256_storeu_pd(&this->data[j * n + k], A);
					}
				}
				for (; k < n; k++) {
					tmp =					this->data[i * n + k];
					this->data[i * n + k] = this->data[j * n + k];
					this->data[j * n + k] = tmp;
				}
			}

			if (pivot != 1.0f) {
				k = 0; tmp = 1/pivot;
				if constexpr (tf32<elem_t> && n > 7) {
					v256s_t A, B = _mm256_set1_ps(tmp);
					for (k = j; k < n; k += 8) {
						A =	_mm256_loadu_ps(&this->data[i * n + k]);
						A = _mm256_mul_ps(A, B);
						_mm256_storeu_ps(&this->data[i * n + k], A);
					}
				} else if constexpr (tf64<elem_t> && n > 3) {
					v256d_t A, B = _mm256_set1_pd(tmp);
					for (k = j; k < n; k += 4) {
						A =	_mm256_loadu_pd(&this->data[i * n + k]);
						A = _mm256_mul_pd(A, B);
						_mm256_storeu_pd(&this->data[i * n + k], A);
					}
				}
				for (; k < n; k++) { this->data[i * n + k] *= tmp; }
			}

			for (j = 0; j < n; j++) {
				if (j == i) { continue; }
				k = 0; tmp = -this->data[j * n + i];
				if constexpr (tf32<elem_t> && n > 7) {
					v256s_t A, B, C = _mm256_set1_ps(tmp);
					for (; k < n; k += 8) {
						A =	_mm256_loadu_ps(&this->data[j * n + k]);
						B =	_mm256_loadu_ps(&this->data[i * n + k]);
						A = _mm256_fmadd_ps(B, C, A);	// B * C + A
						_mm256_storeu_ps(&this->data[j * n + k], A);
					}
				} else if constexpr (tf64<elem_t> && n > 3) {
					v256d_t A, B, C = _mm256_set1_pd(tmp);
					for (; k < n; k += 4) {
						A =	_mm256_loadu_pd(&this->data[j * n + k]);
						B =	_mm256_loadu_pd(&this->data[i * n + k]);
						A = _mm256_fmadd_pd(B, C, A);	// B * C + A
						_mm256_storeu_pd(&this->data[j * n + k], A);
					}
				}
				for (; k < n; k++) {
					this->data[j * m + k] +=  this->data[i * m + k] * tmp;
				}
			}
		}
	}
}


#endif //MATRIX_LIBRARY_H