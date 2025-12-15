//
// Created by marijn on 12/8/25.
//
#include <chrono>
#include <math.h>
#include <stdfloat>

#include "matrix.hpp"


std::chrono::steady_clock::time_point start, end;



void test3x3(void) {
	MAT::sq_matrix<f64_t> mat(3);

	/* 1 2 2
	 * 2 0 2
	 * 2 2 1 */
	mat(0, 0) = 1;
	mat(0, 1) = 2;
	mat(0, 2) = 2;
	mat(1, 0) = 2;
	mat(1, 1) = 0;
	mat(1, 2) = 2;
	mat(2, 0) = 2;
	mat(2, 1) = 2;
	mat(2, 2) = 1;


	mat.print();
	MAT::sq_matrix inverse = mat.inverse();
	//inverse.print();
	mat.RREF();
	mat.print();
}


void test(void) {
	srand(1);
	MAT::sq_matrix<f64_t> mat(1000);
	mat.init_rand();

	//mat.print();
	start = std::chrono::steady_clock::now();
	mat.RREF();
	end = std::chrono::steady_clock::now();
	//mat.print();

	printf("RREF:   %ld us\n", std::chrono::duration_cast<std::chrono::microseconds> (end - start).count());

	srand(1);
	MAT::sq_matrix<f32_t> mat32(1000);
	mat32.init_rand();

	//mat.print();
	start = std::chrono::steady_clock::now();
	mat32.RREF();
	end = std::chrono::steady_clock::now();
	//mat.print();

	printf("RREF32: %ld us\n", std::chrono::duration_cast<std::chrono::microseconds> (end - start).count());


}

int main(void) {
	test3x3();
	test();

	f16_t s = 1.3;

	return 0;
}