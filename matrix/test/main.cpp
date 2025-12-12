//
// Created by marijn on 12/8/25.
//
#include <chrono>

#include "matrix.hpp"



std::chrono::steady_clock::time_point start, end;



void test3x3(void) {
	MAT::sq_matrix mat(3);

	/* 1 2 2
	 * 2 1 2
	 * 2 2 1 */
	mat(0, 0) = 1;
	mat(0, 1) = 2;
	mat(0, 2) = 2;
	mat(1, 0) = 2;
	mat(1, 1) = 1;
	mat(1, 2) = 2;
	mat(2, 0) = 2;
	mat(2, 1) = 2;
	mat(2, 2) = 1;


	//mat.print();
	MAT::sq_matrix inverse = mat.inverse();
	inverse.print();
	//mat.reduce();
	//mat.print();
}


void test_big_rand(void) {
	srand(1);
	MAT::sq_matrix mat(100);
	mat.init_rand();

	MAT::matrix RRWS = mat.get_RRWS();

	mat.print();
	RRWS.print();

	start = std::chrono::steady_clock::now();
	MAT::sq_matrix inv = mat.inverse();
	end = std::chrono::steady_clock::now();

	inv.print();
	RRWS.print();
	printf("calc 100x100 inverse took: %ld ns\n", std::chrono::duration_cast<std::chrono::nanoseconds> (end - start).count());
}


int main(void) {
	//test3x3();
	test_big_rand();  // -> 1000x1000 inverse within 1s

	return 0;
}