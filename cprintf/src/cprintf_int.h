//
// Created by marijn on 11/17/25.
//

#ifndef CPRINTF_CPRINTF_INT_H
#define CPRINTF_CPRINTF_INT_H
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <stdarg.h>

#include "../inc/cprintf.h"



/*!< color_t bit offsets
 * color setting:	bits 0 - 3
 * effect setting:	bits 4 - 6
 * */
#define COLOR_T_CSTART	0
#define COLOR_T_CMASK	(0b1111 << COLOR_T_CSTART)
#define COLOR_T_ESTART	4
#define COLOR_T_EMASK	(0b111 << COLOR_T_ESTART)


// 8 bytes at the start, 4 bytes at the end
#define COLOR_CONST_SIZE	12

// color code range: 30-38
#define COLOR_T_COFFSET		30



__always_inline uint8_t get_color(color_t col)	{ return ((col & COLOR_T_CMASK) >> COLOR_T_CSTART) + COLOR_T_COFFSET; }
__always_inline uint8_t get_effect(color_t col)	{ return ((col & COLOR_T_EMASK) >> COLOR_T_ESTART);}


#endif //CPRINTF_CPRINTF_INT_H
