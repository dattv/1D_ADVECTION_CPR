#include "Include.h"
#include "Element.h"
#include "model.h"

#include "InitialCondition.h"
#include "Method.h"
#include "Output.h"

#include "ArtificialViscousity.h"

#include "Solver.h"
#include "Right_Hand_Side.h"
int main() {

	int order = 5;	int rk_order = 5;
	double time_start = 0.e0;	double time_end = 0.2e0;	double CFL = 5e-2;	double dt = 5e-4;
	
	double type = 0.e0;
	
	NondimensionalNumber<double> noninfo;
	
	CorrectionPolynomial<double> CP_coef;	CP_coef.ComputeCorrectionPolynomial(order);
	
	ComputationalSystem<double> comp;		comp.makeComputationalSystem(order, CP_coef);
	
	MaterialComponents<double> MaterialInfo;
	
	SolverInformation<double> SolverInfo;	SolverInfo.makeSolver(0, rk_order, time_start, time_end, CFL, dt);
	
	Model<double> model;	int nVal	= 2;	double xMin = 0.e0;	double xMax = 1.e0;	int Nele	= 100;
	
	
	/* =========================== Create model with inputs parameter ===========================*/
	model.makeModel(Nele, nVal, xMin, xMax, &comp, &MaterialInfo, &CP_coef, &noninfo, &SolverInfo);
	/* ================================ END CREATION PROCEDURE ==================================*/
	
	
	/* ============================ Compute PI variables for model ==============================*/
	model.Ele->materialInfo->PI = 4.e0*atan(1.e0);
	/* ==========================================================================================*/
	
	
	/* Setup initial condition */
	loop(0, model.nEle, model.Ele, INIT_ADV_1D<decltype(type)>);
	/* ==========================================================================================*/
	
	
	/* Setup boundary element */
	model.BoundaryElement[0] = model.Ele[0];			/* Left boundary */
	model.BoundaryElement[1] = model.Ele[Nele - 1];		/* Right boundary */

	model.Ele[0].left = &model.BoundaryElement[0];
	model.Ele[model.nEle - 1].right = &model.BoundaryElement[1];
	/* ==========================================================================================*/
	
	
	/* Print out initial condition */
	output_1D("initial.tec", model);
	/* ==========================================================================================*/
	
	clock_t begin = clock();
	/* ======================================== work ============================================*/
	//SolverInfo.time_end = 10;
	for (double time = 0; time <= SolverInfo.time_end; time = time + SolverInfo.DT) {

		solver(0, model.nEle, model.Ele, RHS_ADV1D<decltype(type)>);

		int iter = (int)((time - model.Ele[0].solverInfo->time_start) / SolverInfo.DT);
		char fileName[] = "output.tec";
		char sIter[4];
		//if ( iter % 1000 == 0.e0) {
			
			output_1D(fileName, model);
			std::cout << "TIME: " << time << std::endl;
			
		//}

	}
	/* ====================================== End work ===========================================*/
	
	clock_t end = clock();
	std::cout.precision(20);
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "time elapsed: " << elapsed_secs << std::endl;
	
	
	/*============================== Finish program Delete Model ================================*/
	model.deleteModel();
	///* ==========================================================================================*/

	system("pause");
	return 0;
}