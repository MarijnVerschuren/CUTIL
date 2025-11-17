//
// Created by marijn on 11/16/25.
//
#include "cprintf_int.h"



void cprintf(color_t col, const char *__restrict fmt, ...) {
    va_list args;
    va_start(args, fmt);

	uint32_t len = strlen(fmt);
	char* colored_fmt = malloc(strlen(fmt) + COLOR_CONST_SIZE);

	if (col & COLOR_T_EMASK) {
		#ifdef NEGATIVE_NEWLINE_CORRECTION
		uint8_t nnlc = (
			get_effect(col) == get_effect(NEGATIVE)
			&& (fmt[len-1] == '\n')
		);
		if (nnlc) {
			char* temp_fmt = malloc(len);
			memcpy(temp_fmt, fmt, len);
			temp_fmt[len-1] = 0;
			fmt = temp_fmt;
 		}
		#endif
		len = sprintf(
			colored_fmt,
			"\033[%d;%dm%s\033[0m",
			get_color(col),
			get_effect(col),
			fmt
		);
		#ifdef NEGATIVE_NEWLINE_CORRECTION
		if (nnlc) {
			free((void*)fmt);
			colored_fmt[len] = '\n';
			colored_fmt[len+1] = 0;
		}
		#endif
	} else {
		sprintf(
			colored_fmt,
			"\033[0;%dm%s\033[0m",
			get_color(col),
			fmt
		);
	}

	vprintf(colored_fmt, args);
	va_end(args);
	free(colored_fmt);
}