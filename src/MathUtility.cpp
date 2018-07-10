#include "MatUtility.h"

//template <typename T>
//void JacobiGQ(T alpha, T beta, int n, T *x, T *w) {

//	if (n == 0) {
//		x[0] = (alpha - beta) / (alpha + beta + 2.0e0);
//		x[0] = 2.e0;
//	}
//	else {
//		int n1 = n + 1;
//
//		T ireal[n1], h1[n1], jdiag1[n1];
//		for (int i = 0; i < n1; i++) {
//			ireal[i] = (T)i;
//			h1[i] = 2.e0*ireal[i] + alpha + beta;
//			if (h1[i] > 0.e0) {
//				jdiag1[i] = -(alpha*alpha - beta*beta) / (h1[i]*(h1[i] + 2.0e0))
//			}
//			else {
//				jdiga1[i] = 0.e0;
//			}
//		}
//
//		T d1[n], d2[n], d3[n], jdiag2[n];
//		int i1;
//		for (int i = 0; i < n; i++) {
//			d1[i] = 2.e0 / (h1[i] + 2.e0);
//
//			i1 = i + 1;
//			d2[i] = ireal[i1] * (ireal[i1] + alpha + beta)*(ireal[i1] + alpha)*(ireal[i1] + beta);
//
//			d3 = 1.e0 / ((h1[i] + 1.e0)*(h1[i] + 3.e0));
//
//			jdiag2[i] = d1[i] * sqrt(d2[i] * d3[i]);
//		}
//
//
//		int info;
//		T work[max(1, 2 * n)], vect[n1][n1];
//		dsteqr('I', n1, jdiag1, jdiag2, vect, n1, work, info);
//
//
//	}
//	return;
//}