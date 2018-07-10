#pragma once
#ifndef RIGHT_HAND_SIDE_H
#define RIGHT_HAND_SIDE_H

#include "Include.h"
#include "Element.h"
#include "FluxReconstruction.h"

template <typename T>
void RHS_ADV1D(Element<T> *ele) {
	Element<T> *last, *next;

	last = ele->left;
	next = ele->right;

	T phil, phir, ul, ur;
	Vect1D<T> *U;
	int size = ele->U->getSize()[0] - 1;

	if (last != NULL) {
		U = last->U;
		phil = U[0](size);
		ul = U[1](size);
	}
	else {
		U = ele->U;
		phil = U[0](0);
		ul = U[1](0);
	}
	U = ele->U;
	phir = U[0](0);
	ur = U[1](0);

	T NF_l = 0.5e0*(phil*ul + phir*ur - 0.5e0*abs(ul + ur)*(phir - phil));
	T F_l = ur*phir;

	if (next != NULL) {
		U = next->U;
		phir = U[0](0);
		ur = U[1](0);
	}
	else {
		U = ele->U;
		phir = U[0](size);
		ur = U[1](size);
	}
	U = ele->U;
	phil = U[0](size);
	ul = U[1](size);

	T NF_r = 0.5e0*(phil*ul + phir*ur - 0.5e0*(ul + ur)*(phir - phil));
	T F_r = ul*phil;

	U = &ele->U[0];

	CorrectionPolynomial<T> *GB; GB = ele->comSys->CP_coef;
	Vect1D<T> *URHS; URHS = &ele->URHS[0];

	matmul(ele->comSys->Identity, ele->U[1](0), ele->df);

	Fr_Scheme_1D(-1.e0, ele->rx, ele->comSys->Dr, ele->U[0], ele->df, NF_l, NF_r, F_l, F_r, GB->glb, GB->grb, URHS);

	//AuxiliaryU(0, ele);
	//compute_artificialViscosity(0, ele);
	//for (int i = 0; i < URHS->getSize()[0]; i++)
	//URHS->getVec()[i] += ele->q->getVec()[i];

	last	= NULL;
	next	= NULL;
	U		= NULL;
	GB		= NULL;
	URHS	= NULL;

}

template <typename T>
void compute_RHS(int first, int last, Element<T> *ele) {


	Element<T> *lastE, *nextE, *currE;
	int size = ele->U->getSize()[0] - 1;
	Vect1D<T> *U;
	Vect1D<T> DF; DF.makeVec1D(size + 1);
	T phil, phir, ul, ur;
	CorrectionPolynomial<T> *GB;
	Vect1D<T> *URHS;
	T NF_r;
	T F_r;

	for(int iele = first; iele <last; iele ++){
		lastE = ele[iele].left;
		nextE = ele[iele].right;
		currE = &ele[iele];


		U = lastE->U;
		phil = U[0](size);
		ul = U[1](size);

		U = currE->U;
		phir = U[0](0);
		ur = U[1](0);

		T NF_l = 0.5e0*(phil*ul + phir*ur - 0.5e0*abs(ul + ur)*(phir - phil));
		T F_l = ur*phir;

		U = nextE->U;
		phir = U[0](0);
		ur = U[1](0);

		U = currE->U;
		phil = U[0](size);
		ul = U[1](size);

		NF_r = 0.5e0*(phil*ul + phir*ur - 0.5e0*(ul + ur)*(phir - phil));
		F_r = ul*phil;

		U = &currE->U[0];

		GB = currE->comSys->CP_coef;
		URHS = &currE->URHS[0];

		matmul(currE->comSys->Dr, U[0], DF);

		for (int i_order = 0; i_order <= size; i_order++) {
			URHS[0](i_order) = -(DF(i_order) + GB->glb(i_order)*(NF_l - F_l) + GB->grb(i_order)*(NF_r - F_r))*currE->rx(i_order);
		}

	}
	lastE = NULL;
	nextE = NULL;
	U = NULL;
	GB = NULL;
	URHS = NULL;

}
#endif // !RIGHT_HAND_SIDE_H

