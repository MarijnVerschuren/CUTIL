//
// Created by marijn on 11/16/25.
//
#ifndef CPRINTF_CPRINTF_H
#define CPRINTF_CPRINTF_H


typedef enum {
	// colors	[0:3]
	BLACK =		0x00,
	RED =		0x01,
	GREEN =		0x02,
	YELLOW =	0x03,
	BLUE =		0x04,
	PINK =		0x05,
	TEAL =		0x06,
	WHITE =		0x07,
	NONE =		0x08,
	// effects	[4:6]
	BOLD =		0x10,
	FAINT =		0x20,
	CURSIVE =	0x30,
	UNDERLINE =	0x40,
	BLINK =		0x50,
	NEGATIVE =	0x70
} color_t;


#ifdef __cplusplus
color_t operator|(color_t a, color_t b) { return (color_t)(((int)a) | ((int)b)); }
extern "C" {
#endif


void cprintf(color_t col, const char *__restrict fmt, ...);


#ifdef __cplusplus
};
#endif

#endif //CPRINTF_CPRINTF_H
