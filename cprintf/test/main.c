//
// Created by marijn on 11/16/25.
//
#include "cprintf.h"
#include <stdint.h>
#include <stdio.h>


int main() {
	for (uint8_t e = 0; e < 8; e++) {
		for (uint8_t c = 0; c < 8; c++) {
			cprintf(c | (e << 4), "c: %d, e: %d\n", c, e);
		}
		printf("\n");
	}

	return 0;
}