#pragma once
#ifndef ARTIFICIALVISCOUSITY_H
#define ARTIFICIALVISCOUSITY_H

#include "Include.h"
#include "Vect1D.h"
#include "Matrix.h"
#include "MatUtility.h"
#include "FunctionUtility.h"
#include "ComputationalSystem.h"

template<typename T>
void compute_ArtificialViscousityCoef(int id, int first, int last, Element<T> *ele) {

	Element<T> *currE;

	int n = ele->U[0].getSize()[0];
	Vect1D<T> ut;		ut.makeVec1D(n);
	Vect1D<T> u_hat;	u_hat.makeVec1D(n);
	Vect1D<T> uu_hat;	uu_hat.makeVec1D(n);

	Vect1D<T> *cons, *velo;
	T k, so, epso, range1, range2, range3, s, dot_uu, dot_uu_hat, max_val, min_val;

	for (int iele = first; iele < last; iele++) {
		
		currE = &ele[iele];
		cons = &currE->U[id];
		velo = &currE->U[1];

		matmul(currE->comSys->invV, cons[0], ut);
		//currE->comSys->invV.print_out_matrix();
		//cons[0].print_out_vector();
		//ut.print_out_vector();
		ut(n - 1) = 0.e0;

		matmul(currE->comSys->V, ut, u_hat);

		dot_uu = dot_product<T>(cons[0], cons[0]);

		for (int i = 0; i < n; i++)
			uu_hat(i) = cons[0](i) - u_hat(i);


		dot_uu_hat = dot_product(uu_hat, uu_hat);
		//if (abs(dot_uu_hat) <= 1.e-17) dot_uu_hat = 0.e0;
		s = log10(dot_uu_hat / dot_uu);

		if (!isfinite(s) || isnan(s)) s = -6.e0;

		k = 4.5e0;
		so = -4.e0;

		max_val = abs(*std::max_element(velo[0].getVec(), velo[0].getVec() + n));
		min_val = abs(*std::min_element(velo[0].getVec(), velo[0].getVec() + n));

		epso = currE->dx / (cons[0].getSize()[0] - 1)*(max_val > min_val ? max_val : min_val) / 10.e0;

		range1 = abs(s < (so - k));
		range2 = abs((so - k) < s && s < so + k);
		range3 = abs(s >(so + k));

		currE->epsilon[id] = epso*(0.e0*range1 + 0.5e0*(1.e0 + sin(currE->materialInfo->PI*(s - so) / (2.e0*k)))*range2 + range3);
		//std::cout.precision(20);
		//std::cout << currE->ind<< " "<< currE->epsilon[id] << std::endl;
	}
	ut.deleteVec1D();
	u_hat.deleteVec1D();
	uu_hat.deleteVec1D();
}

template <typename T>
void compute_AuxiliaryU(int id, int first, int last, Element<T> *ele) {

	int n1 = ele->U[0].getSize()[0] - 1;
	Element<T> *lastE, *currE, *nextE;
	T k = 0.5e0;
	T Ul, Ur, Vl, Vr;
	Vect1D<T> *U;
	T u_comL;
	T Fl;
	T Fr;
	T u_comR;
	CorrectionPolynomial<T> *GB; 
	Vect1D<T> *tempQ;
	Vect1D<T> DF; DF.makeVec1D(n1 + 1);
	for (int iele = first; iele < last; iele++) {

		lastE = ele[iele].left;
		currE = &ele[iele];
		nextE = ele[iele].right;


		U = &lastE->U[id];
		Ul = U[0](n1);
		Vl = lastE->U[1](n1);


		U = &currE->U[id];
		Ur = U[0](0);
		Vr = currE->U[1](0);

		u_comL = k*Ul*Vl + (1.e0 - k)*Ur*Vr;
		Fl	 = Ur*Vr;


		U = &nextE->U[id];
		Ur = U[0](0);
		Vr = nextE->U[1](0);
		

		Ul = currE->U[id](n1);
		Vl = currE->U[1](n1);
		Fr	 = Ul*Vl;

		u_comR = k*Ul*Vl + (1.e0 - k)*Ur*Vl;

		GB = currE->comSys->CP_coef;
		tempQ = &currE->q[id];
		//Fr_Scheme_1D(currE->epsilon[id], ele->rx, ele->comSys->Dr, ele->U[id], ele->comSys->Identity, u_comL, u_comR, Fl, Fr, GB->glb, GB->grb, tempQ);
		matmul(currE->comSys->Dr, currE->U[id], DF);

		for (int i_order = 0; i_order <= n1; i_order++) {
			tempQ[0](i_order) = currE->epsilon[id] * (DF(i_order) + GB->glb(i_order)*(u_comL - Fl) + GB->grb(i_order)*(u_comR - Fr))*currE->rx(i_order);
		}
		//std::cout << currE->ind;
		//currE->q[id].print_out_vector();
	}
	DF.deleteVec1D();
}

template <typename T>
void compute_artificialViscosity(int id, int first, int last, Element<T> *ele) {

	Element<T> *left, *right, *current;
	T ul, ur, ucom_l, ucom_r, fl, fr;
	T k = 0.5e0;
	int size = ele->q->getSize()[0] - 1;
	Vect1D<T> DF; DF.makeVec1D(size + 1);

	for (int iele = first; iele < last; iele++) {

		left = ele[iele].left;
		current = &ele[iele];
		right = ele[iele].right;

		/* left stage */
		ul = left->q[id](size);

		ur = current->q[id](0);
		ucom_l = k*ul + (1.e0 - k)*ur;
		fl = ur;

		/* right stage */
		ur = right->q[id](0);
		ul = current->q[id](size);
		ucom_r = k*ul + (1.e0 - k)*ur;
		fr = ul;

		matmul(current->comSys->Dr, current->q[id], DF);
		
		for (int i_order = 0; i_order <= size; i_order++) {
			current->Uxx[id](i_order) = (DF(i_order) + current->CP_Coef->glb(i_order)*(ucom_l - fl) + current->CP_Coef->grb(i_order)*(ucom_r - fr))*current->rx(i_order);
		}

		//if (current->ind >=48 ) {
		//	current->comSys->Dr.print_out_matrix();
		//	current->q[id].print_out_vector();
		//	DF.print_out_vector();
		//
		//	std::cout << current->ind;
		//	current->Uxx->print_out_vector();
		//}

	}
	DF.deleteVec1D();
	return;
}

#endif // !ARTIFICIALVISCOUSITY_H

