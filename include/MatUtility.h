#pragma once
#ifndef MATHUTILITY_H
#define MATHUTILITY_H

#include "Include.h"
#include "Vect1D.h"

template <typename T>
void JacobiGQ(T alpha, T beta, int n, T *x, T *w) {

		if (n == 0) {
			x[0] = (alpha - beta) / (alpha + beta + 2.0e0);
			w[0] = 2.e0;
		}
		else {
			int n1 = n + 1;
	
			T *ireal	= new T[n1];
			T *h1		= new T[n1];
			T *jdiag1	=  new T [n1];
			for (int i = 0; i < n1; i++) {
				ireal[i] = (T)i;
				h1[i] = 2.e0*ireal[i] + alpha + beta;
				if (h1[i] > 0.e0) {
					jdiag1[i] = -(alpha*alpha - beta*beta) / (h1[i] * (h1[i] + 2.0e0));
				}
				else {
					jdiag1[i] = 0.e0;
				}
			}
	
			T *d1		= new T[n];
			T *d2		= new T[n];
			T *d3		= new T[n];
			T *jdiag2	= new T [n];

			int i1;
			for (int i = 0; i < n; i++) {
				d1[i] = 2.e0 / (h1[i] + 2.e0);
	
				i1 = i + 1;
				d2[i] = ireal[i1] * (ireal[i1] + alpha + beta)*(ireal[i1] + alpha)*(ireal[i1] + beta);
	
				d3[i] = 1.e0 / ((h1[i] + 1.e0)*(h1[i] + 3.e0));
	
				jdiag2[i] = d1[i] * sqrt(d2[i] * d3[i]);
			}
	
	
			int info;
			T *work = new T [std::max(1, 2 * n)];
			T *vect = new T [n1*n1];
			char *COMPZ = "I";

			dsteqr(COMPZ, &n1, jdiag1, jdiag2, vect, &n1, work, &info);
	
			if (info < 0) {
				std::cout << "Parameter " << n << " in call to dsteqr has illegal value" << std::endl;
			}
			else if (info > 0){
				std::cout << n << " off - diagonal elements have not converged to zero in call to dsteq" << std::endl;
			}
			else {
				for (int i = 0; i < n1; i++) x[i] = jdiag1[i];
				//for (int i = 0; i < n1; i++) std :: cout << vect[i*n1] << " ";
				T lngammaab = lgamma(alpha + beta + 1.e0);
				T lngammaa = lgamma(alpha + 1.e0);
				T lngammab = lgamma(beta + 1.e0);
				for (int i = 0; i < n1; i++) w[i] = vect[i*n1] * vect[i*n1] * pow(2.e0, alpha + beta + 1.e0) / (alpha + beta + 1.e0)*exp(lngammaa + lngammab - lngammaab);
			}
			
			delete[]	ireal	;
			delete[]	h1		;
			delete[]	jdiag1	;
			delete[]	d1		;
			delete[]	d2		;
			delete[]	d3		;
			delete[]	jdiag2  ;
			delete[]	work;
			delete[]	vect	;
		}
		return;
	}

template <typename T>
void JacobiGL(T alpha, T beta, int n, T *JGL) {
	if (n <= 0) {
		std ::cout << "" << std::endl;
		std::cout << "JacobiGL called with n<= 0. Aborting" << std::endl;
	}
	else {
		int n1 = n + 1;
		for (int i = 0; i < n1; i++) JGL[i] = 0.e0;

		if (n == 1) {
			JGL[0] = -1.e0;
			JGL[1] = 1.e0;
		}
		else {
			T *tempJGL = &JGL[1];
			T *w = (T*)mkl_malloc((n - 1) * sizeof(T), 64);// new T[n - 1];
			JacobiGQ(alpha + 1.e0, beta + 1.e0, n - 2, tempJGL, w);
			JGL[0] = - 1.e0;
			JGL[n] = 1.e0;
			mkl_free(w); // delete[] w;
		}
	}
}

template <typename T>
void JacobiP(int sizex, T *x, T alpha, T beta, int n, T *JacobiP) {
	
	int n1 = n + 1;
	T *pl = (T*)mkl_malloc(n1*sizex * sizeof(T), 64);
	

	for (int i = 0; i < n1*sizex; i++) pl[i] = 0.e0;

	T lngammaab = lgamma(alpha + beta + 1.e0);
	T lngammaa = lgamma(alpha + 1.e0);
	T lngammab = lgamma(beta + 1.e0);

	T invSQgamma0 = pow(2.e0, alpha + beta + 1.e0) / (alpha + beta + 1.e0)* exp(lngammaa + lngammab - lngammaab);

	T gamma0 = 1.e0 / sqrt(invSQgamma0);

	for (int i = 0; i < sizex; i++) pl[i] = gamma0;
	
	if (n == 0) {
		for (int i = 0; i < sizex; i++) JacobiP[i] = pl[i];
	}
	else {
		T gamma1 = 0.5e0*sqrt((alpha + beta + 3.e0) / ((alpha + 1.e0)*(beta + 1.e0)))*gamma0;
		T fac1 = (alpha + beta + 2.e0);
		T fac2 = (alpha - beta);
		for (int i = 0; i < sizex; i++) pl[sizex + i] = gamma1 *(fac1*x[i] + fac2);

		if (n == 1) {
			for (int i = 0; i < sizex; i++) JacobiP[i] = pl[sizex + i];
		}
		else {
			T aold = 2.e0 / (2.e0 + alpha + beta)*sqrt((1.e0 + alpha)*(1.e0 + beta) / (3.e0 + alpha + beta));

			for (int i = 0; i < n - 1; i++) {
				T ireal = (T)i + 1.e0;
				T irealp1 = ireal + 1.e0;
				T h1 = 2.e0*ireal + alpha + beta;
				T anew = 2.e0 / (h1 + 2.e0)*sqrt(irealp1*(irealp1 + alpha + beta) *(irealp1 + alpha)*(irealp1 + beta) / (h1 + 1.e0) / (h1 + 3.e0));
				T bnew = -(pow(alpha, 2.e0) - pow(beta, 2.e0)) / (h1*(h1 + 2.e0));
				
				for (int j = 0; j < sizex; j++) pl[(i + 2)*sizex + j] = 1.e0 / anew *(-aold*pl[i*sizex + j] + (x[j] - bnew)*pl[(i + 1)*sizex + j]);
				aold = anew;
			}
			for (int j = 0; j < sizex; j++) JacobiP[j] = pl[n*sizex + j];
			//std::cout << " " << std::endl;
			//for (int j = 0; j < n1; j++) {
			//	for (int k = 0; k < sizex; k++) std::cout << pl[j*sizex + k] << " ";
			//	std::cout << " " << std::endl;
			//}
		}
	}

	mkl_free(pl);
	return;
}

template <typename T>
void GradJacobiP(int sizex, T *x, T alpha, T beta, int n, T *GJP) {

	if (n == 0) {
		for (int i = 0; i < sizex; i++) GJP[i] = 0.e0;
	}
	else {
		JacobiP<T>(sizex, x, alpha + 1.e0, beta + 1.e0, n - 1, GJP);
		for (int i = 0; i < sizex; i++) GJP[i] *= sqrt(n*(n + alpha + beta + 1.e0));
	}
	return;
}

template <typename T>
void Vandemonde1D(int sizex, int n, T *x, T *V) {
	int n1 = n + 1;
	T *tempV = (T*)mkl_malloc(sizex * sizeof(T), 64);
	for (int j = 0; j < n1; j++) {
		JacobiP<T>(sizex, x, 0.e0, 0.e0, j, tempV);
		for (int k = 0; k < sizex; k++) V[k*n1 + j] = tempV[k];
	}
	mkl_free(tempV);
}

template <typename T>
void GradVandermonde1D(int sizex, int n, T *x, T* DV) {

	int n1 = n + 1;
	T *tempDV = (T*)mkl_malloc(sizex * sizeof(T), 64);
	for (int j = 0; j < n1; j++) {
		GradJacobiP<T>(sizex, x, 0.e0, 0.e0, j, tempDV);
		for (int k = 0; k < sizex; k++) DV[k*n1 + j] = tempDV[k];
	}

	mkl_free(tempDV);
	return;
}

template <typename T>
void Dmatrix1D(int n, T *x, T*V, T*Dr) {

	int n1 = n + 1;

	T *vt	= (T*)mkl_malloc((n1)*(n1) * sizeof(T), 64);
	T *vrt  = (T*)mkl_malloc((n1)*(n1) * sizeof(T), 64);
	T *vr	= (T*)mkl_malloc((n1)*(n1) * sizeof(T), 64);
	
	transpose<T>(n1, V, vt);

	GradVandermonde1D<T>(n1, n, x, vr);

	transpose<T>(n1, vr, vrt);

	int *ipiv = (int*)mkl_malloc(n1 * sizeof(int), 64);
	int info;

	info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n1, n1, vt, n1, ipiv, vrt, n1);
	//dgesv(&n1, &n1, vt, &n1, ipiv, vrt, &n1, &info);

	transpose(n1, vrt, Dr);

	mkl_free(vt);
	mkl_free(vrt);
	mkl_free(vr);
	mkl_free(ipiv);
	return;
}

template <typename T>
void transpose(int n, T *m1, T *m2) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			m2[j*n + i] = m1[i*n + j];
		}
	}
	return;
}

/* Auxiliary routine: printing a matrix */
template <typename T>
void print_matrix(char* desc, int m, int n, T* a) {
	int i, j;
	printf("\n %s\n", desc);
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) printf(" %10.4f", a[i*n + j]);
		printf("\n");
	}
}

template <typename T>
void inverseMat(int size,T * A,T* invA) {
	int *ipiv = (int*)mkl_malloc(size * sizeof(int), 64);
	T *work = (T*)mkl_malloc(size * sizeof(T), 64);
	int info;

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) invA[i*size + j] = A[i*size + j];
	}

	//LAPACKE_dgetrf(LAPACK_ROW_MAJOR, size, size, invA, size, ipiv, info);
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, size, size, invA, size, ipiv);
	//LAPACKE_dgetri(LAPACK_ROW_MAJOR, size, invA, size, ipiv, work, size, info);
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, size, invA, size, ipiv);
	mkl_free(ipiv);
	mkl_free(work);
}

template <typename T>
void Jacobian(int n,T* vx, T *r, T *dr, T *x, T *rx) {

	int n1 = n + 1;
	for (int i = 0; i < n1; i++) {
		x[i] = vx[0] + 0.5e0*(1.e0 + r[i])*(vx[1] - vx[0]);
	}

	T xr;
	for (int i = 0; i < n1; i++) {
		xr = 0.e0;
		for (int j = 0; j < n1; j++) {
			xr += dr[i*n1 + j] * x[j];
		}
		rx[i] = 1.e0 / xr;
	}

}

template <class T>
T dot_product(Vect1D<T> &v1, Vect1D<T> &v2) {
	T sum = 0.e0;
	//if (std::is_same(T, double)) {
		sum = cblas_ddot(v1.getSize()[0], v1.getVec(), 1, v2.getVec(), 1);
	//}
	//else if (std::is_same(T, float)) {
	//	sum = cblas_sdot(v1.getSize()[0], v1.getVec(), 1, v2.getVec(), 1);
	//}
	return sum;
}

template <class T>
void subtract(int n, T *v1, T *v2, T* res) {
	for (int i = 0; i < n; i++) res[i] = v1[i] - v2[i];
	return;
}
#endif // !MATHUTILITY_H

