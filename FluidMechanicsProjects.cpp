// FluidMechanicsProjects.cpp : Defines the entry point for the console application.
//
// VR: this is the main interface use for calling all the different apps (unit testing)

#include "stdafx.h"
#include <iostream> 
#include "REE\RadialEquilibriumEquationSolver.h"
#include "FileParser.h"

int _tmain(int argc, _TCHAR* argv[])
{
	std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
	std::cout << "Reading input file" << std::endl;
	FileParser fileParser;
	std::string myPath("C:\\Users\\Valentin\\Documents\\Coding\\FluidMechanicsProjects\\REE\\input\\problem2.dat");
	std::map<std::string, std::string> myMap;
	myMap = fileParser.readFile(myPath);
	if (myMap.empty())
	{
		std::cout << "Problem encountered when reading input file" << std::endl;
	}

	std::cout << "Entering main app" << std::endl;
	REESolver reeSolver;
	reeSolver.setData(myMap);
	reeSolver.runREESolver();
	//reeSolver.run1DmeanlineAnalysis();
	return 0;
}

