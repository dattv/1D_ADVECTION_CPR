#pragma once
#ifndef ELEMENT_H
#define ELEMENT_H

#include "Include.h"
#include "MaterialComponents.h"
#include "ComputationalSystem.h"
#include "CorrectionPolynomial.h"
#include "NondimensionalNumber.h"
#include "SolverInfomation.h"
#include "Title.h"

template <class T>
class Element : public Title
{
public:
	Element();
	~Element();

	void MakeElement(int order, int nval, int ind, T xl, T xr, ComputationalSystem<T> *comSys,
														MaterialComponents<T> *materialInfo,
														CorrectionPolynomial<T> *CP_Coef,
														NondimensionalNumber<T> *dimensionalInfo,
														SolverInformation<T> *solverInfo);
	void deleteElement();

	int order;
	int ind;
	int nval;
	
	Element<T> *left  = NULL;
	Element<T> *right = NULL;

	T vx[2];
	T xc;
	T dx;
	Vect1D<T> xk;

	ComputationalSystem<T> *comSys;

	Vect1D<T> rx; /* jacobian */

	T localTime;

	Vect1D<T> *U;	/* conservative  */

	Vect1D<T> *URHS;

	Vect1D<T> *ResU;

	T *epsilon;

	Vect1D<T> *q; /*auxilirary*/
	Vect1D<T> *Uxx;

	Matrix<T> df;

	T *AV;

	MaterialComponents<T> *materialInfo;

	CorrectionPolynomial<T> *CP_Coef;

	NondimensionalNumber<T> *dimensionalInfo;

	SolverInformation<T> *solverInfo;

};

template <class T>
Element<T>::Element()
{
	this->name = "Element";
}


template <class T>
Element<T>::~Element()
{
	this->name = "";
}

template <class T>
void Element<T>::MakeElement(int order, int nval, int ind, T xl, T xr, ComputationalSystem<T> *comSys,
																MaterialComponents<T> *materialInfo,
																CorrectionPolynomial<T> *CP_Coef,
																NondimensionalNumber<T> *dimensionalInfo,
																SolverInformation<T> *solverInfo) {
	this->name = "Element";
	this->order = order;
	this->ind = ind;
	this->nval = nval;

	this->vx[0] = xl;
	this->vx[1] = xr;

	this->dx = (xr - xl);
	this->xc = 0.5e0*this->dx;

	this->localTime = 0.e0;

	this->U		= new Vect1D<T> [nval];
	this->URHS	= new Vect1D<T> [nval];
	this->ResU = new Vect1D<T>[nval];

	this->epsilon = new T [nval];
	this->q = new Vect1D<T>[nval];
	this->Uxx = new Vect1D<T>[nval];

	Vect1D<T> *temPTR;
	for (int i = 0; i < nval; i++) {
		temPTR = &this->U[i];
		temPTR->makeVec1D(order + 1);

		temPTR = &this->URHS[i];
		temPTR->makeVec1D(order + 1);

		temPTR = &this->ResU[i];
		temPTR->makeVec1D(order + 1);

		temPTR = &this->q[i];
		temPTR->makeVec1D(order + 1);

		temPTR = &this->Uxx[i];
		temPTR->makeVec1D(order + 1);
	}
	temPTR = NULL;

	this->AV = (T*)mkl_malloc(nval * sizeof(T), 64);

	this->comSys = comSys;
	this->materialInfo = materialInfo;
	this->CP_Coef = CP_Coef;
	this->dimensionalInfo = dimensionalInfo;
	this->solverInfo = solverInfo;

	ComputationalSystem<T> *tempComSys = this->comSys;
	//MaterialComponents<T> *tempMaterialInfo = this->materialInfo;

	this->xk.makeVec1D(order + 1);
	this->rx.makeVec1D(order + 1);

	Jacobian(this->order, this->vx, tempComSys->r.getVec(), tempComSys->Dr.getMatrix(), this->xk.getVec(), this->rx.getVec());

	this->df.makeMatrix(order + 1, order + 1);
}

template <class T>
void Element<T> ::deleteElement() {
	this->name = "";

	this->order = 0;
	this->ind = 0;
	this->nval = 0;

	this->vx[0] = 0;
	this->vx[1] = 0;

	this->xc = 0;
	this->dx = 0;

	this->localTime = 0.e0;

	Vect1D<T> *tempPTR;
	for (int i = 0; i < nval; i++) {
		tempPTR = &this->U[i];
		tempPTR->deleteVec1D();

		tempPTR = &this->URHS[i];
		tempPTR->deleteVec1D();

		tempPTR = &this->ResU[i];
		tempPTR->deleteVec1D();

		tempPTR = &this->q[i];
		tempPTR->deleteVec1D();

		tempPTR = &this->Uxx[i];
		tempPTR->deleteVec1D();
	}

	delete [] this->U;
	delete [] this->URHS;
	delete [] this->ResU;

	delete[] this->epsilon;

	delete[] this->q;
	
	delete[] this->Uxx;

	mkl_free(this->AV );

	this->comSys = NULL;
	this->materialInfo = NULL;
	this->CP_Coef = NULL;
	this->dimensionalInfo = NULL;
	this->solverInfo = NULL;

	this->xk.deleteVec1D();
	this->rx.deleteVec1D();

	this->df.deleteMatrix();
}
#endif // !ELEMENT_H

