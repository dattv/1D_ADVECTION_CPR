#pragma once
#ifndef RIEMANNINFORMATION_H
#define RIEMANNINFORMATION_H

#include "Include.h"
#include "Title.h"

template<class T>
class RiemannInformation : public Title
{
public:
	RiemannInformation();
	~RiemannInformation();


	int RiemannType = 0;
	T alpha = 0.5e0;
};

template<class T>
RiemannInformation<T>::RiemannInformation()
{
	this->name = "RiemannInfo";
}

template<class T>
RiemannInformation<T>::~RiemannInformation()
{
	this->name = "RiemannInfo";
}
#endif // !RIEMANNINFORMATION_H

