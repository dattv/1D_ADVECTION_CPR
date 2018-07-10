#pragma once
#ifndef SOLVERINFORMATION_H
#define SOLVERINFORMATION_H

#include "Include.h"
#include "Title.h"
#include "Vect1D.h"

template <class T>
class SolverInformation : public Title
{
public:
	SolverInformation();
	~SolverInformation();

	void makeSolver(int solver_type, int nStep, T t_start, T t_end, T CFL, T DT);

	bool lSolve = 1;
	int	 solver_type; /* 0 ruku; 
					  1 explicit second order, 
					  2 implicit, 
					  3 semi-imlicit, 
					  4 pressure-based semi-implicit */
	int  nstep;
	Vect1D<T> rka;
	Vect1D<T> rkb;
	Vect1D<T> rkc;

	T alpha = 0.5e0;
	T time_start = 0.e0;
	T time_end = 0.e0;
	T globalTime = 0.e0;
	T CFL = 0.05e0;
	T DT = 0.e0;

};

template <class T>
SolverInformation<T>::SolverInformation()
{
}

template <class T>
SolverInformation<T>::~SolverInformation()
{
}

template <class T>
void SolverInformation<T>::makeSolver(int solver_type, int nStep, T t_start, T t_end, T CFL, T DT) {

	this->CFL = CFL;
	this->DT = DT;
	this->time_start = t_start;
	this->time_end = t_end;

	if (solver_type == 0) {
		if (nStep != 5) {
			std::cout << "error";
			return;
		}
		this->nstep = nStep;
		this->name = "EX_Runge_Kutta";
		this->rka.makeVec1D(nStep);
		this->rkb.makeVec1D(nStep);
		this->rkc.makeVec1D(nStep);

		this->nstep = nStep;

		this->rka(0) = 0.0e0;
		this->rka(1) = -567301805773.0e0 / 1357537059087.0e0;
		this->rka(2) = -2404267990393.0e0 / 2016746695238.0e0;
		this->rka(3) = -3550918686646.0e0 / 2091501179385.0e0;
		this->rka(4) = -1275806237668.0e0 / 842570457699.0e0;

		this->rkb(0) = 1432997174477.0e0 / 9575080441755.0e0;
		this->rkb(1) = 5161836677717.0e0 / 13612068292357.0e0;
		this->rkb(2) = 1720146321549.0e0 / 2090206949498.0e0;
		this->rkb(3) = 3134564353537.0e0 / 4481467310338.0e0;
		this->rkb(4) = 2277821191437.0e0 / 14882151754819.0e0;

		this->rkc(0) = 0.0e0;
		this->rkc(1) = 1432997174477.0e0 / 9575080441755.0e0;
		this->rkc(2) = 2526269341429.0e0 / 6820363962896.0e0;
		this->rkc(3) = 2006345519317.0e0 / 3224310063776.0e0;
		this->rkc(4) = 2802321613138.0e0 / 2924317926251.0e0;

	}
	else if (solver_type == 1) {
		this->name = "EX_Second_Order";
	}
	else if (solver_type == 2) {
		this->name = "Implicit";
	}
	else if (solver_type == 3) {
		this->name = "Semi_implicit";
	}
	else if (solver_type == 4) {
		this->name = "Pressure_Base_Semi_implicit";
	}
		
}
#endif // !SOLVERINFORMATION_H

