//
// Created by marijn on 12/8/25.
//
#include "matrix.hpp"


int main(void) {
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


	mat.print();
	mat.reduce();
	mat.print();

	return 0;
}