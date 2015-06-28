// FluidMechanicsProjects.cpp : Defines the entry point for the console application.
//
// VR: this is the main interface use for calling all the different apps (unit testing)

#include "stdafx.h"
#include <iostream> 
#include "REE\RadialEquilibriumEquationSolver.h"

int _tmain(int argc, _TCHAR* argv[])
{
	std::cout << "Entering main app\n" << std::endl;
	REESolver reeSolver;
	reeSolver.runREESolver();
	reeSolver.run1DmeanlineAnalysis();
	return 0;
}

