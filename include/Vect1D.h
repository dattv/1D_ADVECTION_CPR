#pragma once
#ifndef VECT1D_H
#define VECT1D_H

#include "Include.h"
#include "Title.h"

template <class T>
class Vect1D : public Title
{
public:
	Vect1D();
	~Vect1D();

	void makeVec1D(int size);
	void deleteVec1D();

	int *getSize();
	T *getVec();


	Vect1D<T> operator+ (Vect1D<T> &v1);
	Vect1D<T> operator+ (T val);

	Vect1D<T> operator- (Vect1D<T> &v1);
	Vect1D<T> operator- (T val);

	Vect1D<T> operator* (Vect1D<T> &v1);
	Vect1D<T> operator* (T val);

	Vect1D<T> operator/ (Vect1D<T> &v1);
	Vect1D<T> operator/ (T val);

	Vect1D<T> operator= (Vect1D<T> &v1);
	Vect1D<T> operator= (T val);

	T &operator () (int ind);

	Vect1D<T> subtract(Vect1D<T> &v1, Vect1D<T> &v2);

	void print_out_vector();

private:
	int size = 0;
	T *array;
};

template <class T>
Vect1D<T>::Vect1D()
{
}

template <class T>
Vect1D<T>::~Vect1D()
{
}

template <class T>
void Vect1D<T> :: makeVec1D(int size) {
	this->name = "Vec1D";
	this->size = size;
	this->array = (T*)mkl_malloc(size * sizeof(T), 64);
	//this->array = new T[size];
	for (int i = 0; i < size; i++)
		this->array[i] = 0.e0;
}

template <class T>
void Vect1D<T> ::deleteVec1D() {
	this->name = "";
	this->size = 0;
	mkl_free(this->array);
	//delete[] this->array;
}


template <class T>
int *Vect1D<T>::getSize() {
	return &this->size;
}

template <class T>
T *Vect1D<T>::getVec() {
	return this->array;
}

template <class T>
Vect1D<T> Vect1D<T>::operator+ (Vect1D<T> &v1) {
	if (v1.getSize()[0] != this->getSize()[0]) {
		std::cout << v1.getSize()[0] << this->getSize()[0];
	}
	else {
		for (int i = 0; i < this->getSize()[0]; i++) {
			this->array[i] += v1.getVec()[i];
		}
	}
	return *this;
}

template <class T>
Vect1D<T> Vect1D<T>::operator- (Vect1D<T> &v1) {
	if (v1.getSize()[0] != this->getSize()[0]) {
		std::cout << v1.getSize()[0] << this->getSize()[0];
	}
	else {
		for (int i = 0; i < this->getSize()[0]; i++) {
			this->array[i] -= v1.getVec()[i];
		}
	}
	return *this;
}

template <class T>
Vect1D<T> Vect1D<T>::operator + (T val) {
	for (int i = 0; i < this->getSize()[0]; i++) {
		this[i] = val;
	}
}

template <class T>
Vect1D<T> Vect1D<T>::operator- (T val) {
	for (int i = 0; i < this->getSize()[0]; i++) {
		this[i] -= val;
	}
}

template <class T>
Vect1D<T> Vect1D<T>::operator* (Vect1D<T> &v1) {
	for (int i = 0; i < this->getSize()[0]; i++) {
		this[i] *= v1[i];
	}
}

template <class T>
Vect1D<T> Vect1D<T>::operator/ (Vect1D<T> &v1) {
	for (int i = 0; i < this->getSize()[0]; i++) {
		this[i] /= v1[i];
	}
}

template <class T>
Vect1D<T> Vect1D<T>::operator/ (T val) {
	for (int i = 0; i < this->getSize()[0]; i++) {
		this[i] /= val;
	}
}

template <class T>
Vect1D<T> Vect1D<T>::operator* (T val) {
	for (int i = 0; i < this->getSize()[0]; i++) {
		this->getVec()[i] *= val;
	}
	return *this;
}

template <class T>
Vect1D<T> Vect1D<T>::operator= (Vect1D<T> &v1) {
	for (int i = 0; i < this->getSize()[0]; i++) {
		this->getVec()[i] = v1.getVec()[i];
		//this[i] = v1[i];
	}
	return *this;
}

template <class T>
Vect1D<T> Vect1D<T>::operator= (T val) {
	for (int i = 0; i < this->getSize()[0]; i++) {
		this[i] = val;
	}
}

template <class T>
T &Vect1D<T>::operator () (int ind) {
	return this->array[ind];
}

template <class T>
void Vect1D<T>::print_out_vector() {
	std::cout << " " << std::endl;
	std::cout.precision(20);
	for (int i = 0; i < this->getSize()[0]; i++) {
		std::cout << this->array[i] << "    ";
	}
	std::cout << " " << std::endl;
}

template <class T>
Vect1D<T> Vect1D<T>::subtract(Vect1D<T> &v1, Vect1D<T> &v2) {
	Vect1D<T> tempVec;
	tempVec.makeVec1D(v1.getVec()[0]);
	for (int i = 0; i < v1.getVec()[0]; i++) tempVec(i) = v1(i) - v2(i);
	return  tempVec;
}
#endif // !VECT1D_H

