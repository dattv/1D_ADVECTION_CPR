#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include "Include.h"
#include "Element.h"
#include "Method.h"
#include "RungekuttaSolver.h"

template<typename T>
void intialResidual_single(Element<T> *ele) {

	for (int i = 0; i < ele->ResU->getSize()[0]; i++)
		ele->ResU[0](i) = 0.e0;
}

template <typename T, typename Functional>
int solver(int first, int last, Element<T> *ele, Functional RHS) {

	T templateType = 0.e0;
	/* initial residual */
	loop(first, last, ele, intialResidual_single<decltype(templateType)>);

	if (ele->solverInfo->solver_type == 0) { /* Explicit Runge Kutta */
		EX_RungeKutta(first, last, ele);
	}
	else if (ele->solverInfo->solver_type == 1) {/* Explicit second order*/

	}
	else if (ele->solverInfo->solver_type == 2) {/* Implicit */
	
	}
	else {
		std::cout << "Not supported this solver yet" << std::endl;
		return 1;
	}


	return 0;
}

#endif // !SOLVER_H

