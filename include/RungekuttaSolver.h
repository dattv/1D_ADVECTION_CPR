#pragma once
#ifndef RUNGEKUTTASOLVER_H
#define RUNGEKUTTASOLVER_H

#include "Include.h"
#include "Element.h"
#include "Right_Hand_Side.h"
#include "Method.h"

template <typename T>
void EX_RungeKutta(int first, int last, Element<T> *ele) {
	Element<T> *tempEle;

	for (int i = 0; i < ele->solverInfo->nstep; i++) { /* loop over runge kutta */

		ele->localTime = ele->solverInfo->globalTime + ele->solverInfo->rkc(i) * ele->solverInfo->DT;

		/* compute Right Hanside */
		//loop(first, last, ele, RHS);
		compute_RHS(first, last, ele);

		compute_ArtificialViscousityCoef(0, first, last, ele);
		compute_AuxiliaryU(0, first, last, ele);

		///* compute AuxiliaryU */
		//for (int iele = first; iele < last; iele++) {
		//	tempEle = &ele[iele];
		//	AuxiliaryU(0, tempEle);
		//}
		

		/* compute Artificial viscosity */
		//for (int iele = first; iele < last; iele++) {
		//	tempEle = &ele[iele];
		//	compute_artificialViscosity(0, tempEle);
		//	//std::cout << "====================================";
		//	//tempEle->q[0].print_out_vector();
		//	//tempEle->Uxx[0].print_out_vector();
		//}
		//
		compute_artificialViscosity(0, first, last, ele);
		/* adding Artificial viscosity to RHS */
		for (int iele = first; iele < last; iele++) {
			tempEle = &ele[iele];
			for (int j = 0; j < tempEle->URHS->getSize()[0]; j++) {
				tempEle->URHS[0](j) += tempEle->Uxx[0](j);
			}
			//std::cout << tempEle->ind;
			//tempEle->URHS[0].print_out_vector();
			//std::cout << "=====================================";
		}

		for (int j = first; j < last; j++) { /* loop from fist element to last element*/
			tempEle = &ele[j];
			Vect1D<T> *ResU, *Urhs, *U;
			for (int k = 0; k < tempEle->nval; k++) { /* loop in number of variables */

				ResU = &tempEle->ResU[k];
				Urhs = &tempEle->URHS[k];
				U = &tempEle->U[k];

				for (int i_order = 0; i_order < ResU->getSize()[0]; i_order++) { /* loop inside node of each element */

					ResU[0](i_order) = ele->solverInfo->rka(i) * ResU[0](i_order)
						- ele->solverInfo->DT *Urhs[0](i_order);

					U[0](i_order) -= ele->solverInfo->rkb(i) * ResU[0](i_order);
					//if (U[0](i_order) < 0.e0) U[0](i_order) = 0.e0;
					//if (U[0](i_order) > 1.e0) U[0](i_order) = 1.e0;
				}
			}
			//if (j == 50) {
			//std::cout << tempEle->ind;
			//tempEle->ResU[0].print_out_vector();
			//tempEle->U[0].print_out_vector();
			//std::cout << " ";
			//}

		}
	}

	for (int i = first; i < last; i++) {
		tempEle = &ele[i];
		tempEle->solverInfo->globalTime += tempEle->solverInfo->DT;
		tempEle->localTime = tempEle->solverInfo->globalTime;
	}
}


#endif // !RUNGEKUTTASOLVER_H

