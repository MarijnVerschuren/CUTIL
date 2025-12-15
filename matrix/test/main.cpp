//
// Created by marijn on 12/8/25.
//
#include <chrono>
#include <math.h>

#include "matrix.hpp"



std::chrono::steady_clock::time_point start, end;



void test3x3(void) {
	MAT::sq_matrix mat(3);

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
	mat.REF();
	mat.print();
}


void test_4x4a(void) {
	MAT::matrix amat(4, 4);
	amat(0, 0) = -9;				//amat(0, 0) = 1;
	amat(0, 1) = -36;				//amat(0, 1) = 0;
	amat(0, 2) = 27;				//amat(0, 2) = 1;
	amat(0, 3) = -36;				//amat(0, 3) = 0;
	//amat(0, 4) = 0;				//amat(0, 4) = 0;
	amat(1, 0) = -6;				//amat(1, 0) = 2;
	amat(1, 1) = -23;				//amat(1, 1) = 1;
	amat(1, 2) = 22;				//amat(1, 2) = 1;
	amat(1, 3) = -28;				//amat(1, 3) = 1;
	//amat(1, 4) = 0;				//amat(1, 4) = 0;
	amat(2, 0) = 5;				//amat(2, 0) = 2;
	amat(2, 1) = 17;				//amat(2, 1) = 2;
	amat(2, 2) = -27;				//amat(2, 2) = 1;
	amat(2, 3) = 32;				//amat(2, 3) = 1;
	//amat(2, 4) = 0;				//amat(2, 4) = sqrt(3);
	amat(3, 0) = 2;				//amat(3, 0) = 0;
	amat(3, 1) = 7;				//amat(3, 1) = 2;
	amat(3, 2) = -10;				//amat(3, 2) = 0;
	amat(3, 3) = 12;				//amat(3, 3) = 1;
	//amat(3, 4) = 0;				//amat(3, 4) = sqrt(3) - 1;
	amat.print();
	amat.REF();
	amat.print();
}


void test_big_rand(void) {
	srand(1);
	MAT::sq_matrix mat(1000);
	mat.init_rand();

	MAT::matrix RRWS = mat.get_RRWS();

	//mat.print();
	//RRWS.print();

	start = std::chrono::steady_clock::now();
	MAT::sq_matrix inv = mat.inverse();
	end = std::chrono::steady_clock::now();

	inv.print();

	//RRWS.print();
	printf("calc 100x100 inverse took: %ld ns\n", std::chrono::duration_cast<std::chrono::nanoseconds> (end - start).count());
}


void test(void) {
	srand(1);
	MAT::sq_matrix mat(1000);
	mat.init_rand();

	start = std::chrono::steady_clock::now();
	mat.REF_C();
	end = std::chrono::steady_clock::now();

	printf("REF_C: %ld ns\n", std::chrono::duration_cast<std::chrono::nanoseconds> (end - start).count());

	mat.init_rand();
	start = std::chrono::steady_clock::now();
	mat.REF();
	end = std::chrono::steady_clock::now();

	printf("REF: %ld ns\n", std::chrono::duration_cast<std::chrono::nanoseconds> (end - start).count());


}

int main(void) {
	test3x3();
	//test_4x4a();	// used for systems q5!
	//test_big_rand();  // -> 1000x1000 inverse within 1s
	test();

	return 0;
}