#include <iostream>
#include <vector>
#include <fstream>

#include "fem_matrix.h"
#include "parameter.h"
#include "output.h"
#include "structure_solver.h"
#include "eigen_solver.h"
#include "time_measure.h"
#include "init.h"
#include "debug.h"

void executeStaticAnalysis(){
	// read input files
	readParameterDataFile();
	readNodeDataFile();
	readElemDataFile();
	readBoundaryDataFile();
    // set initial state and memory
    initField();
    initAdjMatrix();
	// output initial state
	outputParameterDataFile();
	outputDispVtkFile(0);
	// calculate boundary shape function
	calcBoundaryShapeFunction();
	// coloring Elements for Element Matrix Parallelization
	coloringElements();
	
#ifdef MEASURE
	std::ofstream ofs(sim_prm.time_output_filename, std::ios::app);
#endif
	// initialize sparse matrixs
	initSparseMatrix();
	if(sim_prm.dim == 2){
		calcElementMatrix2Dquad();
	}else if(sim_prm.dim == 3){
		calcElementMatrix3Dtetra();
	}else{
		std::cout<<"Dimension Error: "<<sim_prm.dim<<std::endl;
		exit(1);
	}
	// set sparse matrix
	setSparseMatrix();
	// solve KU=F
	solveLinearEquation2D();
	// output
	outputDispVtkFile(1);

#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}