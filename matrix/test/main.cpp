//
// Created by marijn on 12/8/25.
//
#include <chrono>
#include <math.h>
#include <stdfloat>

#include "complex.hpp"
#include "matrix.hpp"


std::chrono::steady_clock::time_point start, end;



// void test(void) {
// 	srand(1);
// 	MAT::sq_matrix<f64_t> mat(1000);
// 	mat.init_rand();
//
// 	//mat.print();
// 	start = std::chrono::steady_clock::now();
// 	mat.RREF();
// 	end = std::chrono::steady_clock::now();
// 	//mat.print();
// 	printf("RREF:   %ld us\n", std::chrono::duration_cast<std::chrono::microseconds> (end - start).count());
//
// 	mat.init_rand();
// 	start = std::chrono::steady_clock::now();
// 	mat *= 1.23456789;
// 	end = std::chrono::steady_clock::now();
//
// 	printf("mul assign:   %ld us\n", std::chrono::duration_cast<std::chrono::microseconds> (end - start).count());
//
// 	srand(1);
// 	MAT::sq_matrix<f32_t> mat32(1000);
// 	mat32.init_rand();
//
// 	//mat.print();
// 	start = std::chrono::steady_clock::now();
// 	mat32.RREF();
// 	end = std::chrono::steady_clock::now();
// 	//mat.print();
//
// 	printf("RREF32: %ld us\n", std::chrono::duration_cast<std::chrono::microseconds> (end - start).count());
//
//
// }


void test_math(void) {
	MAT::matrix<f64_t,4,4> mat_a;
	MAT::matrix<f64_t,4,4> mat_b;
	mat_a.init_rand();
	mat_b.init_rand();
	mat_a.print();
	mat_b.print();

	MAT::matrix<f64_t,4,4> mat_c = mat_a + mat_b;
	mat_c.print();

	mat_c *= 8;
	mat_c.print();
	mat_c.REF();
	mat_c.print();

	mat_b += mat_a * -5;
	mat_b.print();
	mat_b.RREF();
	mat_b.print();

	MAT::matrix<f64_t,4,4> mat_d = mat_a * mat_b;
}

void test_mat_mul() {
	MAT::matrix<f64_t,4,5> mat_a;	mat_a.init_rand();
	MAT::matrix<f64_t,5,3> mat_b;	mat_b.init_rand();
	MAT::matrix<f64_t,4,3> mat_c =	mat_a * mat_b;
	mat_a.print();
	mat_b.print();
	mat_c.print();

}


void test_timing() {
	MAT::matrix<f64_t, 1000, 1000> mat_a;	mat_a.init_rand();
	MAT::matrix<f64_t, 1000, 1000> mat_b;	mat_b.init_rand();

	start = std::chrono::steady_clock::now();
	MAT::matrix<f64_t, 1000, 1000> mat_c = mat_a * mat_b;
	end = std::chrono::steady_clock::now();
	printf("1000x1000 mul took: %llu ns\n", std::chrono::duration_cast<std::chrono::microseconds> (end - start).count());
}

void test_C() {
	MAT::complex<f64_t> c1(3, 6);
	MAT::complex<f64_t> c2(-2, -1.42);
	c1.print();
	c2.print();
	MAT::complex<f64_t> c3 = c1 + c2;
	c3.print();
	c3 *= 6;
	c3.print();
	((MAT::matrix<f64_t, 2, 2>)c3).print();
}


int main(void) {
	srand(1);
	//test_math();
	//test_mat_mul();

	//test_timing();
	test_C();

//	f16_t s = 1.3;

	return 0;
}