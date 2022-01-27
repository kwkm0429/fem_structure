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

void solve(){
	int t, idx=1;

	// read input files
	readParameterDataFile();
	readNodeDataFile();
	readElemDataFile();
	readBoundaryDataFile();
    std::cout<<"Succeeded Initialization\n"<<std::endl;

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
	
	for(t=1; t<=sim_prm.max_time_step; t++){
#ifdef MEASURE
		std::ofstream ofs(sim_prm.time_output_filename, std::ios::app);
		ofs<<"-------------------- TIME STEP "<<t<<" --------------------"<<std::endl;
#endif
		std::cout<<"Time step: "<<t<<std::endl;

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
		// solve static analysis
		solveStaticAnalysis();

		if(t%sim_prm.output_interval == 0){
			outputDispVtkFile(idx);
			idx++;
		}
		std::cout<<std::endl;
	}
}