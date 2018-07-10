#pragma once
#ifndef OUTPUT_H
#define OUTPUT_H

#include "Include.h"
#include "model.h"

template <class T>
void output_1D(std::string fileName, Model<T> model) {
	std::ofstream myFile;

	myFile.open(fileName);
	Element<T> *E;
	if (myFile.is_open()) {
		for (int i = 0; i < model.nEle; i++) {
			E = &model.Ele[i];
			for (int j = 0; j < E->order; j++) {
				//std::cout << E->xk(j);
				myFile << E->xk(j);
				for (int k = 0; k < E->nval; k++) {
					myFile << "   " << E->U[k](j);
				}
				myFile << " \n";
			}

		}
	}
}
#endif // !OUTPUT_H

