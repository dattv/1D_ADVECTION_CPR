#pragma once
#ifndef METHOD_H
#define METHOD_H

#include "Include.h"
#include "model.h"


template <typename T, typename functional>
void implement(int ind, Element<T> *ele, functional f) {

	f(&ele[ind]);
	return;
}

template <class T, typename functional>
void  loop(int first, int last, Element<T> *ele, functional function) {

	
	for (int i = first; i < last; i++) {
		//implement(i, ele, function);
		implement(i, ele, function);
	}

}

#endif // !METHOD_H

