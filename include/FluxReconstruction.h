#pragma once
#ifndef FLUXRECONSTRUCTION_H
#define FLUXRECONSTRUCTION_H

#include "Include.h"
#include "Vect1D.h"
#include "Matrix.h"

template <typename T>
void Fr_Scheme_1D(T epsilon, Vect1D<T> J, Matrix<T> Dr, Vect1D<T> U, Matrix<T> A, T FN_l, T FN_r, T F_l, T F_r, Vect1D<T> G_lb, Vect1D<T> G_rb, Vect1D<T> *res) {
	Vect1D<T> temp;
	temp.makeVec1D(U.getSize()[0]);

	matmul(Dr, U, res[0]);
	matmul(res[0], A, temp);
	//res[0].print_out_vector();
	for (int i = 0; i < res->getSize()[0]; i++) {
		res[0](i) =  epsilon*(temp(i)+ G_lb(i)*(FN_l - F_l) + G_rb(i)*(FN_r - F_r))*J(i);

	}
	temp.deleteVec1D();
	return;
}
#endif // !FLUXRECONSTRUCTION_H

