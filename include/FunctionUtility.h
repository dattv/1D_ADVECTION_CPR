#pragma once
#ifndef FUNCTIONUTILITY_H
#define FUNCTIONUTILITY_H

#include"Include.h"
#include "Vect1D.h"
#include "Matrix.h"

template <typename T>
void matmul(Matrix<T> m1, T m2, Matrix<T> &m3) {
	int m_m1 = m1.getMsize()[0];
	int n_m1 = m1.getNsize()[0];

	//Matrix<T> tempMat;	tempMat.makeMatrix(m_m1, n_m1);

	for (int i = 0; i < m_m1; i++) {
		for (int j = 0; j < n_m1; j++)
			m3(i,j) = m1(i,j)*m2;
	}
	return;
}

template <typename T>
void matmul(Matrix<T> m1, Vect1D<T> m2, Vect1D<T> &v) {
	int m_m1 = m1.getMsize()[0];
	int n_m1 = m1.getNsize()[0];

	int m_m2 = m2.getSize()[0];
	int n_m2 = 1;

	//Vect1D<T> v;	v.makeVec1D(m_m1);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_m1, n_m2, n_m1, 1.e0, m1.getMatrix(), n_m1, m2.getVec(), n_m2, 0.e0, v.getVec(), n_m2);
	//v.deleteVec1D();
	//T temp;
	//T *Arr1 = m1.getMatrix();
	//T *Arr2 = m2.getVec();
	//T *res = v.getVec();
	//for (int i = 0; i < m_m1; i++) {
	//	temp = 0.e0;
	//	for (int j = 0; j < n_m1; j++)
	//		temp += Arr1[i*n_m1 + j] * Arr2[j];
	//	res[i] = temp;
	//}
	//Arr1 = NULL;
	//Arr2 = NULL;
	//res = NULL;
	return; 
}

template <typename T>
void matmul(Vect1D<T> m1, Matrix<T> m2, Vect1D<T> &v) {
	int m_m1 = 1;
	int n_m1 = m1.getSize()[0];

	int m_m2 = m2.getMsize()[0];
	int n_m2 = m2.getNsize()[0];

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_m1, n_m2, n_m1, 1.e0, m1.getVec(), n_m1, m2.getMatrix(), n_m2, 0.e0, v.getVec(), n_m2);
	//v.print_out_vector();
	//T temp;
	//T *Arr1 = m1.getVec();
	//T *Arr2 = m2.getMatrix();
	//T *res = v.getVec();
	//for (int i = 0; i < n_m2; i++) {
	//	temp = 0.e0;
	//	for (int j = 0; j < n_m1; j++) {
	//		temp += Arr1[j] * Arr2[i*n_m1 + j];
	//	}
	//	res[i] = temp;
	//}
	return ;
}
#endif // !FUNCTIONUTILITY_H

