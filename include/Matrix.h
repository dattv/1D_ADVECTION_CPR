#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include "Include.h"
#include "Vect1D.h"

template <class T>
class Matrix : public Title
{
public:
	Matrix();
	~Matrix();

	int *getMsize();
	int *getNsize();

	T & operator () (int row, int col);
	void makeMatrix(int nSize, int mSize);
	void makeIdentity(int nSize, int mSize);

	T *getMatrix();

	Matrix<T> operator + (Matrix<T> &m1);
	Matrix<T> operator + (T val);

	Matrix<T> operator - (Matrix<T> &m1);
	Matrix<T> operator - (T val);

	Matrix<T> operator * (Matrix<T> &m1);
	Matrix<T> operator * (T val);

	Matrix<T> operator & (Matrix<T> &m1);


	Matrix<T> operator / (Matrix<T> &m1);
	Matrix<T> operator / (T val);


	Matrix<T> operator = (Matrix<T> &m1);
	Matrix<T> operator = (T val);

	Matrix<T> matmul(Matrix<T> m1, Matrix<T> m2);
	Matrix<T> matmul(Matrix<T> m1, Vect1D<T> m2);

	Matrix<T> transpose();
	void deleteMatrix();

	void print_out_matrix();

private:
	
	int mSize;
	int nSize;
	Vect1D<T> Array;
};

template <class T>
Matrix<T>::Matrix()
{
}

template <class T>
Matrix<T>::~Matrix()
{
}

template <class T>
void Matrix<T>::makeMatrix(int mSize, int nSize) {
	this->name = "Matrix";
	this->nSize = nSize;
	this->mSize = mSize;

	this->Array.makeVec1D(nSize*mSize);
}

template <class T>
void Matrix<T>::deleteMatrix() {
	this->name = " ";
	this->nSize = 0;
	this->mSize = 0;
	this->Array.deleteVec1D();
}

template <class T>
T & Matrix<T>::operator () (int row, int col) {
	return this->Array(this->mSize*(row) + col);
}

template <class T>
int *Matrix<T>::getMsize() {
	return &this->mSize;
}

template <class T>
int *Matrix<T>::getNsize() {
	return &this->nSize;
}

template <class T>
Matrix<T> Matrix<T>::operator + (Matrix<T> &m1) {
	for (int i = 0; i < m1.getMsize()[0]; i++) {
		for (int j = 0; j < m1.getNsize()[0]; j++) {
			this->Array[i*this->getNsize()[0] + j] += m1(i, j);
		}
	}
}

template <class T>
Matrix<T> Matrix<T>::operator + (T val) {
	for (int i = 0; i < this->getMsize()[0]; i++) {
		for (int j = 0; j < this->getNsize()[0]; j++) {
			this->Array[i*this->getNsize()[0] + j] += val;
		}
	}
}

template <class T>
Matrix<T> Matrix<T>::operator - (Matrix<T> &m1) {
	for (int i = 0; i< this->getMsize()[0]; i++) {
		for (int j = 0; j < this->getNsize()[0]; j++) {
			this->Array[i*this->getNsize()[0] + j] -= m1(i, j);
		}
	}
}

template <class T>
Matrix<T> Matrix<T>::operator - (T val) {
	for (int i = 0; i< this->getMsize()[0]; i++) {
		for (int j = 0; j < this->getNsize()[0]; j++) {
			this->Array[i*this->getNsize()[0] + j] -= val;
		}
	}
}

template <class T>
Matrix<T> Matrix<T>::operator * (Matrix<T> &m1) {
	for (int i = 0; i< this->getMsize()[0]; i++) {
		for (int j = 0; j < this->getNsize()[0]; j++) {
			this->Array(i*this->getNsize()[0] + j) *= m1(i, j);
		}
	}
	return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator * (T val) {
	for (int i = 0; i< this->getMsize()[0]; i++) {
		for (int j = 0; j < this->getNsize()[0]; j++) {
			this->Array[i*this->getNsize()[0] + j] *= val;
		}
	}
}

template <class T>
Matrix<T> Matrix<T>::operator / (Matrix<T> &m1) {
	for (int i = 0; i< this->getMsize()[0]; i++) {
		for (int j = 0; j < this->getNsize()[0]; j++) {
			this->Array[i*this->getNsize()[0] + j] /= m1(i, j);
		}
	}
}

template <class T>
Matrix<T> Matrix<T>::operator / (T val) {
	for (int i = 0; i< this->getMsize()[0]; i++) {
		for (int j = 0; j < this->getNsize()[0]; j++) {
			this->Array[i*this->getNsize()+j] /= val;
		}
	}
}

template <class T>
Matrix<T> Matrix<T>::operator = (Matrix<T> &m1) {
	if (this->getMsize()[0] != m1.getMsize()[0] && this->getMsize()[0] != m1.getNsize()[0]) {
		this->makeMatrix(m1.getMsize()[0], m1.getNsize()[0]);
	}
	for (int i = 0; i< this->getMsize()[0]; i++) {
		for (int j = 0; j < this->getNsize()[0]; j++) {
			this->Array(i*this->getNsize()[0] + j) = m1(i, j);
		}
	}
	return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator = (T val) {
	for (int i = 0; i< this->getMsize()[0]; i++) {
		for (int j = 0; j < this->getNsize()[0]; j++) {
			this->Array(i*getNsize()[0] + j) = val;
		}
	}
	return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator & (Matrix<T> &m1) {
	Matrix<T> temp;
	T tempval;
	temp = *this;

	for (int i = 0; i < this->getMsize()[0]; i++) {
		for (int j = 0; j < m1.getNsize()[0]; j++) {
			tempval = 0.e0;
			for (int k = 0; k < this->getNsize()[0]; k++) {
				tempval += temp(i, k) *m1(k, j);
			}
			this->Array(i*this->getNsize()[0] + j) = tempval;
		}

	}
	return *this;
}

template <class T>
void Matrix<T>::print_out_matrix() {
	std::cout.precision(20);
	std::cout << " " << std::endl;
	for (int i = 0; i < this->getMsize()[0]; i++) {
		for (int j = 0; j < this->getNsize()[0]; j++) {
			std::cout << this->Array(i*this->getNsize()[0] + j) << "   ";
		}
		std::cout << " " << std::endl;
	}
	std::cout << " " << std::endl;
}

template <class T>
T *Matrix<T>::getMatrix() {
	return this->Array.getVec();
}


template <class T>
Matrix<T> Matrix<T>::matmul(Matrix<T> m1, Matrix<T> m2) {

	int m_m1 = *m1.getMsize();
	int n_m1 = *m1.getNsize();

	int m_m2 = *m2.getMsize();
	int n_m2 = *m2.getNsize();

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_m1, n_m2, n_m1, 1.e0, m1.getMatrix(), n_m1, m2.getMatrix(), n_m2, 0.e0, this->getMatrix(), n_m2);
	return *this;
}

template <class T>
Matrix<T> Matrix<T>::matmul(Matrix<T> m1, Vect1D<T> v) {

	int m_m1 = m1.getMsize()[0];
	int n_m1 = m1.getNsize()[0];

	int m_m2 = v.getSize()[0];
	int n_m2 = 1;
	//m1.print_out_matrix();
	//v.print_out_vector();

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_m1, n_m2, n_m1, 1.e0, m1.getMatrix(), n_m1, v.getVec(), n_m2, 0.e0, this->getMatrix(), n_m2);
	return *this;
}

template <class T>
Matrix<T> Matrix<T>::transpose() {
	int m = this->getMsize()[0];
	int n = this->getNsize()[0];

	Matrix<T> temp;
	temp.makeMatrix(m, n);
	temp = *this;

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			this->getMatrix()[i*n + j] = temp.getMatrix()[j*m + i];
		}
	}
	this->getMsize()[0] = n;
	this->getNsize()[0] = m;
	temp.deleteMatrix();
	return *this;
}

template <class T>
void Matrix<T>::makeIdentity(int nSize, int mSize) {
	this->makeMatrix(nSize, mSize);
	for (int i = 0; i < nSize; i++) {
		for (int j = 0; j < mSize; j++) {
			if (i == j)
				this->Array.getVec()[i*mSize + j] = 1.e0;
		}
	}
}
#endif // !MATRIX_H

