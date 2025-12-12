//
// Created by marijn on 2/12/23.
//

#ifndef STD_LIB_MATH_H
#define STD_LIB_MATH_H

#include "defs.h"


// arithmetic.asm
extern f64_t fabs(f64_t x);
extern f32_t fabs32(f32_t x);

// trig.asm
extern f64_t sin(f64_t x);
extern f64_t cos(f64_t x);
extern f64_t tan(f64_t x);
extern f64_t csc(f64_t x);
extern f64_t sec(f64_t x);
extern f64_t cot(f64_t x);

// exp.asm
extern f64_t sqrt(f64_t x);
extern f64_t ln(f64_t x);
extern f64_t log2(f64_t x);
extern f64_t log10(f64_t x);
extern f64_t exp(f64_t x);
extern f64_t pow(f64_t x, f64_t y);


#endif //STD_LIB_MATH_H
