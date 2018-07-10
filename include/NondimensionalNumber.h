#pragma once
#ifndef NONDIMENSIONALNUMBER_H
#define NONDIMENSIONALNUMBER_H

#include "Include.h"
#include "Title.h"

template <class T>
class NondimensionalNumber : public Title
{
public:
	NondimensionalNumber();
	~NondimensionalNumber();


	T cInfty	 = 1.e0;
	T rhoInfty	 = 1.e0;
	T TInfty	 = 1.e0;
	T muInfty	 = 1.e0;
	T cpInfty	 = 1.e0;
	T kappaInfty = 1.e0;
	T lInfty	 = 1.e0;
	T pInfty = 1.e0;
	T EInfty = 1.e0;

};

template <class T>
NondimensionalNumber<T>::NondimensionalNumber()
{
	this->name = "NondimensionalNumbers";
}


template <class T>
NondimensionalNumber<T>::~NondimensionalNumber()
{
	this->name = "";
}
#endif // !NONDIMENSIONALNUMBER_H

