//
// Created by marijn on 12/8/25.
//
#include <chrono>
#include <math.h>
#include <stdfloat>

#include "complex.hpp"
#include "matrix.hpp"


std::chrono::steady_clock::time_point start, end;


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
	mat_d.print();
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
	MAT::matrix<f32_t, 1000, 1000> mat_a;	mat_a.init_rand();
	MAT::matrix<f32_t, 1000, 1000> mat_b;	mat_b.init_rand();

	start = std::chrono::steady_clock::now();
	MAT::matrix<f32_t, 1000, 1000> mat_c = mat_a * mat_b;
	end = std::chrono::steady_clock::now();
	printf("1000x1000 mul took: %llu ns\n", std::chrono::duration_cast<std::chrono::microseconds> (end - start).count());
}

void test_C() {
	MAT::complex<f64_t> c1(2, 2);
	MAT::complex<f64_t> c2(4, 1);
	MAT::complex<f64_t> c3 = c1 * c2;
	c1.print();
	c2.print();
	c3.print();
	MAT::sq_matrix<f64_t, 2> m3 = c1.mat();
	m3.print();
	MAT::sq_matrix<f64_t, 2> m3i = m3.inverse();
	m3i.print();
	MAT::complex<f64_t> c4(m3i);
	c4.print();
}


int main(void) {
	srand(1);

	test_math();
	test_timing();
	test_C();

	return 0;
}