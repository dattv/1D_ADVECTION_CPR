#pragma once
#ifndef COMPUTATIONALSYSTEM_H
#define COMPUTATIONALSYSTEM_H

#include "Include.h"
#include "Title.h"
#include "Vect1D.h"
#include "Matrix.h"
#include "MatUtility.h"
#include "CorrectionPolynomial.h"

template <class T>
class ComputationalSystem	:	public Title
{
public:
	ComputationalSystem();
	~ComputationalSystem();

	void makeComputationalSystem(int order, CorrectionPolynomial<T> &CP_coef);

	int order;
	Vect1D<T> r;
	Matrix<T> V;
	Matrix<T> invV;
	Matrix<T> Dr;

	CorrectionPolynomial<T> *CP_coef;

	Matrix<T> Identity;
};

template <class T>
ComputationalSystem<T>::ComputationalSystem()
{
}

template <class T>
ComputationalSystem<T>::~ComputationalSystem()
{
}

template <class T>
void ComputationalSystem<T>::makeComputationalSystem(int order, CorrectionPolynomial<T> &CP_coef) {
	this->order = order;
	
	this->r.makeVec1D(order + 1);
	this->V.makeMatrix(order + 1, order + 1);
	this->invV.makeMatrix(order + 1, order + 1);
	this->Dr.makeMatrix(order + 1, order + 1);

	JacobiGL(0.e0, 0.e0, order, this->r.getVec());

	Vandemonde1D(this->r.getSize()[0], order, this->r.getVec(), this->V.getMatrix());
	
	Dmatrix1D(order, this->r.getVec(), this->V.getMatrix(), this->Dr.getMatrix());
	
	inverseMat(order + 1, this->V.getMatrix(), this->invV.getMatrix());
	
	//this->invV.print_out_matrix();
	
	this->CP_coef = &CP_coef;
	
	this->Identity.makeIdentity(order + 1, order + 1);
	return;
}
#endif // !COMPUTATIONALSYSTEM_H

