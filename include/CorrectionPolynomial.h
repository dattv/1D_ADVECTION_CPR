#pragma once
#ifndef CORRECTIONPOLYNOMIAL_H
#define CORRECTIONPOLYNOMIAL_H

#include "Include.h"
#include "Title.h"
#include "Vect1D.h"

template <class T>
class CorrectionPolynomial : public Title
{
public:
	CorrectionPolynomial();
	~CorrectionPolynomial();

	void ComputeCorrectionPolynomial(int order);

	int order;
	Vect1D<T> glb;
	Vect1D<T> grb;
};

template <class T>
CorrectionPolynomial<T>::CorrectionPolynomial()
{
}


template <class T>
CorrectionPolynomial<T>::~CorrectionPolynomial()
{
}

template <class T>
void CorrectionPolynomial<T>::ComputeCorrectionPolynomial(int order) {
	this->name = "CorrectionPolynomial";
	this->order = order;
	this->glb.makeVec1D(order + 1);
	this->grb.makeVec1D(order + 1);

	if (order == 0) {
		this->glb(0) = -0.5e0;
	}
	else if (order == 1) {
		this->glb(0) = -2.e0; 
			this->glb(1) = 1.e0	;
	}
	else if (order == 2) {
		this->glb(0) = -4.5e0;
		this->glb(1) = 0.75e0;
		this->glb(2) = -1.5e0;
	}
	else if (order == 3) {
		this->glb(0) = -8.e0;
		this->glb(1) = 0.894427190999916e0;
		this->glb(2) = -0.894427190999916e0;
		this->glb(3) = 2.000000000000000e0;
	}
	else if (order == 4) {
		this->glb(0) = -12.500000000000000e0;
		this->glb(1) = 1.071428571428571e0;
		this->glb(2) = -0.937500000000000e0;
		this->glb(3) = 1.071428571428571e0;
		this->glb(4) = -2.500000000000000e0;
	}
	else if (order == 5) {
		this->glb(0) = -18.000000000000000e0;
		this->glb(1) = 1.259090802393860e0;
		this->glb(2) = -1.039883175166253e0;
		this->glb(3) = 1.039883175166253e0;
		this->glb(4) = -1.259090802393862e0;
		this->glb(5) = 3.000000000000000e0;
	}
	else if (order == 6) {

		this->glb(0) = -24.500000000000000e0;
		this->glb(1) = 1.451626611323438e0;
		this->glb(2) = -1.162370412976332e0;
		this->glb(3) = 1.093750000000000e0;
		this->glb(4) = -1.162370412976332e0;
		this->glb(5) = 1.451626611323440e0;
		this->glb(6) = -3.500000000000000e0;
	}
	else if (order == 7) {
		this->glb(0) = -32.000000000000000e0;
		this->glb(1) = 1.646813245455385e0;
		this->glb(2) = -1.294272340806527e0;
		this->glb(3) = 1.177038562354133e0;
		this->glb(4) = -1.177038562354133e0;
		this->glb(5) = 1.294272340806526e0;
		this->glb(6) = -1.646813245455381e0;
		this->glb(7) = 4.000000000000000e0;
	}
	else if (order == 8) {
		this->glb(0) = -40.500000000000000e0;
		this->glb(1) = 1.843607008227080e0;
		this->glb(2) = -1.431394905143639e0;
		this->glb(3) = 1.274249435378094e0;
		this->glb(4) = -1.230468750000000e0;
		this->glb(5) = 1.274249435378094e0;
		this->glb(6) = -1.431394905143638e0;
		this->glb(7) = 1.843607008227083e0;
		this->glb(8) = -4.500000000000002e0;
	}
	else if (order == 9) {
		this->glb(0) = -50.000000000000000e0;
		this->glb(1) = 2.041450781795803e0;
		this->glb(2) = -1.571734950051846e0;
		this->glb(3) = 1.379242344610341e0;
		this->glb(4) = -1.302362052293982e0;
		this->glb(5) = 1.302362052293982e0;
		this->glb(6) = -1.379242344610340e0;
		this->glb(7) = 1.571734950051850e0;
		this->glb(8) = -2.041450781795814e0;
		this->glb(9) = 4.999999999999996e0;
	}
	else if (order == 10) {
		this->glb(0) = -60.500000000000000e0;
		this->glb(1) = 2.240019259474176e0;
		this->glb(2) = -1.714207384864337e0;
		this->glb(3) = 1.489064083963912e0;
		this->glb(4) = -1.384624051179488e0;
		this->glb(5) = 1.353515625000000e0;
		this->glb(6) = -1.384624051179487e0;
		this->glb(7) = 1.489064083963914e0;
		this->glb(8) = -1.714207384864336e0;
		this->glb(9) = 2.240019259474186e0;
		this->glb(10) = -5.500000000000004e0;
	}
	else if (order == 11) {
		this->glb(0) = -72.000000000000000e0;
		this->glb(1) = 2.439109405612857e0;
		this->glb(2) = -1.858170261712427e0;
		this->glb(3) = 1.602105877713413e0;
		this->glb(4) = -1.473343869993111e0;
		this->glb(5) = 1.417653739264046e0;
		this->glb(6) = -1.417653739264047e0;
		this->glb(7) = 1.473343869993112e0;
		this->glb(8) = -1.602105877713410e0;
		this->glb(9) = 1.858170261712421e0;
		this->glb(10) = -2.439109405612830e0;
		this->glb(11) = 6.000000000000000e0;
	}
	else {
		return;
	}

	for (int i = 0; i < order + 1; i++) {
		this->grb(i) = -this->glb(order - i);
	}
}
#endif // !CORRECTIONPOLYNOMIAL_H

