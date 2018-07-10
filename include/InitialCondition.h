#pragma once

#ifndef INITIALCONDITION_H
#define INITIALCONDITION_H

#include "Include.h"
#include "Element.h"


template <typename T>
void INIT_ADV_1D(Element<T> *ele) {

	Vect1D<T> *PHI = &ele->U[0];

	for (int i = 0; i < ele->order + 1; i++) {
		PHI[1](i) = 1.e0;
		if (ele->xk(i) <= 0.5e0)
			PHI[0](i) = 1.e0;
		else
			PHI[0](i) = 0.e0;
	}
}


#endif // !INITIALCONDITION_H


