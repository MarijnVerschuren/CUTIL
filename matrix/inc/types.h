//
// Created by marijn on 12/11/25.
//
#ifndef TYPES_H
#define TYPES_H
#include <stdint.h>
#include <immintrin.h>

typedef uint8_t		u8_t;
typedef uint16_t	u16_t;
typedef uint32_t	u32_t;
typedef uint64_t	u64_t;
typedef int8_t		i8_t;
typedef int16_t		i16_t;
typedef int32_t		i32_t;
typedef int64_t		i64_t;

typedef __m128		v128s_t;	// [f32, f32, f32, f32]
typedef __m128d		v128d_t;	// [f64, f64]
typedef __m128d		v128i_t;	// 16x u8_t ... 2x u64_t
typedef __m256		v256s_t;	// [f32, f32, f32, f32, f32, f32, f32, f32]
typedef __m256d		v256d_t;	// [f64, f64, f64, f64]
typedef __m256d		v256i_t;	// 32x u8_t ... 4x u64_t

typedef double		f64_t;
typedef float		f32_t;
typedef _Float16	f16_t;

#endif //TYPES_H
