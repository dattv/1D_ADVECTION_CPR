#pragma once
#ifndef MATERIALCOMPONENTS_H
#define MATERIALCOMPONENTS_H

#include "Include.h"
#include "Title.h"

template <class T>
class MaterialComponents : public Title
{
public:
	MaterialComponents();
	~MaterialComponents();

	T  gamma = 1.4e0;
	T  B = 1.e0;
	T  rho = 1.3e0;
	T  R_gas_constant = 8.31457e0;
	T  mu = 18.27e0;
	T  lambda = 1.512041288e0;
	T  C_sutherland_constant = 120.e0;
	T  temperature = 291.15e0;

	T PI = 3.141592653589793;
};

template <class T>
MaterialComponents<T>::MaterialComponents()
{
	this->name = "MaterialComponents";
}

template <class T>
MaterialComponents<T>::~MaterialComponents()
{
	this->name = "";
}
#endif // !MATERIALCOMPONENTS_H

