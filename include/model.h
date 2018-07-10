#pragma once
#ifndef MODEL_H
#define MODEL_H

#include "Include.h"
#include "Element.h"
#include "Title.h"

template<class T>
class Model : public Title
{
public:
	Model();
	~Model();

	void makeModel(int nele, int nval, T lB, T rB, ComputationalSystem<T> *comSys,
											MaterialComponents<T> *materialInfo,
											CorrectionPolynomial<T> *CP_Coef,
											NondimensionalNumber<T> *dimensionalInfo,
											SolverInformation<T> *solverInfo);

	void deleteModel();

	Element<T> *Ele;

	Element<T> *BoundaryElement;
	int nEle;
	int nBoundaryElement;

};

template<class T>
Model<T>::Model()
{
	this->name = "Model1D";
}


template<class T>
Model<T>::~Model()
{
	this->name = "";
}

template <class T>
void Model<T>::makeModel(int nEle, int nval, T lB, T rB, ComputationalSystem<T> *comSys,
															MaterialComponents<T> *materialInfo,
															CorrectionPolynomial<T> *CP_Coef,
															NondimensionalNumber<T> *dimensionalInfo,
															SolverInformation<T> *solverInfo) {
	
	this->Ele = new Element<T> [nEle];
	this->nEle = nEle;
	T dx = (rB - lB) / nEle;
	T xl, xr;
	int order = comSys->order;
	for (int i = 0; i < nEle; i++) {
		xl = i*dx;
		xr = xl + dx;
		
		Ele[i].MakeElement(order, nval, i, xl, xr, comSys, materialInfo, CP_Coef, dimensionalInfo, solverInfo);
		if (i < nEle - 1) Ele[i].right = &Ele[i + 1];
		if (i > 0) Ele[i].left = &Ele[i - 1];

	}

	/* Recalculate the DT*/
	this->Ele->solverInfo->DT = this->Ele->solverInfo->CFL*dx / 1.e0;


	// create boundary element
	this->nBoundaryElement = 2;
	this->BoundaryElement = new Element<T>[nBoundaryElement];

}

template <class T>
void  Model<T>::deleteModel() {

	for (int i = 0; i < this->nEle; i++) {
		Ele[i].deleteElement();
	}
	delete[] this->Ele;
}
//
//template <class T>
//void Model<T>::implement(int eleInd, void(function)(Element<T> *)) {
//
//	Element<T> *elePTR = &this->Ele[eleInd];
//	function(elePTR);
//
//	return;
//}
//
//template <class T>
//void Model<T>::loop(int firstEleInd, int lastEleInd, void(function)(Element<T> *)) {
//
//	if (firstEleInd >= 0 && lastEleInd < this->nEle) {
//		for (int i = firstEleInd; i <= lastEleInd; i++) this->implement(i, function(this->Ele));
//	}
//	return;
//}
//
//template <class T>
//void Model<T>::loop_full(void(function)(Element<T> *)) {
//	
//	loop(0, this->nEle - 1, function(this->Ele));
//	return;
//}
#endif // !MODEL_H

