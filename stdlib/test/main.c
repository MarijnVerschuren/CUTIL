//
// Created by marijn on 2/12/23.
//
#include <stdio.h>
#include <malloc.h>

#include "cpu.h"
#include "math.h"
#include "string.h"
#include "vector.h"


void basic_trig(void) {
	f64_t x;
	for (uint8_t i = 0;; i++) {		// loop that roughly divides pi into 255 pieces
		x = (f64_t)i / 0x51; printf("%d:\n", i);
		printf("sin(%f)\t=> %.20f\n", x, sin(x));
		printf("cos(%f)\t=> %.20f\n", x, cos(x));
		printf("tan(%f)\t=> %.20f\n", x, tan(x));
		printf("csc(%f)\t=> %.20f\n", x, csc(x));
		printf("sec(%f)\t=> %.20f\n", x, sec(x));
		printf("cot(%f)\t=> %.20f\n", x, cot(x));
		if (i == 0xff) { break; }  // check at the end to include 0xff
	}
}
void basic_exp(void) {
	f64_t x;
	for (uint8_t i = 0;; i++) {		// loop that roughly divides pi into 255 pieces
		x = (f64_t)i / 0xff; printf("%d:\n", i);
		printf("sqrt(%d)\t\t\t\t=> %.20f\n", i, sqrt(i));
		printf("ln(%f)\t\t\t=> %.20f\n", x, ln(x));
		printf("log2(%f)\t\t\t=> %.20f\n", x, log2(x));
		printf("log10(%f)\t\t\t=> %.20f\n", x, log10(x));
		printf("exp(%f)\t\t\t=> %.20f\n", x, exp(x));
		printf("pow(%f, %f)\t=> %.20f\n", x, x, pow(x, x));
		if (i == 0xff) { break; }  // check at the end to include 0xff
	}
}
void cpuid_test(void) {
	uint32_t a = 0, b = 0, c = 0, d = 0;
	for (uint32_t leaf = 0; leaf < 0x17; leaf++) {
		cpuid(leaf, &a, &b, &c, &d);
		printf("%-20lx => A: %-20luB: %-20luC: %-20luD: %lu\n", leaf, a, b, c, d);
	}	printf("\n");
	for (uint32_t leaf = 0x80000000; leaf < 0x80000022; leaf++) {
		cpuid(leaf, &a, &b, &c, &d);
		printf("%-20lx => A: %-20luB: %-20luC: %-20luD: %lu\n", leaf, a, b, c, d);
	}	printf("\n");
	cpuid(AMD_EASTER_EGG, &a, &b, &c, &d);
	printf("%-20s => A: %-20luB: %-20luC: %-20luD: %lu\n", "AMD_EASTER_EGG", a, b, c, d);
	char* brand_string = calloc(49, 1);
	cpuid(BRAND_STRING0, (uint32_t*)brand_string, (uint32_t*)(brand_string + 4), (uint32_t*)(brand_string + 8), (uint32_t*)(brand_string + 12));
	cpuid(BRAND_STRING1, (uint32_t*)(brand_string + 16), (uint32_t*)(brand_string + 20), (uint32_t*)(brand_string + 24), (uint32_t*)(brand_string + 28));
	cpuid(BRAND_STRING2, (uint32_t*)(brand_string + 32), (uint32_t*)(brand_string + 36), (uint32_t*)(brand_string + 40), (uint32_t*)(brand_string + 44));
	printf("%-20s => %s\n", "BRAND_STRING", brand_string);
	char* vendor_id_string = calloc(13, 1);
	cpuid(VENDOR_ID, (uint32_t*)vendor_id_string, (uint32_t*)vendor_id_string, (uint32_t*)(vendor_id_string + 8), (uint32_t*)(vendor_id_string + 4));
	printf("%-20s => %s\n", "VENDOR_ID", vendor_id_string);
}
void string_test(void) {
	char test[30] = "hello sssssssssssssssssss!\n";
	printf("%lld\n", strlen("hello world!\n"));
	printf("%lld\n", strlen(test));
	strcpy(test, "hello world!\n");
	strcat(test, "hello computer!\n");
	printf("%lld -> %s\n", strlen(test), test);

	char text[] = "a b c d";
	printf("%s\n", strtok(text, " "));
	printf("%s\n", strtok(NULL, NULL));
	printf("%s\n", strtok(NULL, NULL));
	printf("%s\n", strtok(NULL, NULL));
}

int main(int argc, char** argv) {
	//basic_trig();
	//basic_exp();
	cpuid_test();
	//string_test();
	//test();

	return 0;
}
