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
//#include "debug.h"

void updatePosition(Str& str){
	int i;
	for(int i;i<str.num_nodes;i++){
		str.x[i] += str.disp_x[i];
		str.y[i] += str.disp_y[i];
		str.z[i] += str.disp_z[i];
	}
}

void readInputFiles(Sim& sim, Str& str){
	// read input files
	readParameterDataFile(sim, str);
	readNodeDataFile(sim, str);
	readElemDataFile(sim, str);
	readBoundaryDataFile(sim, str);
	for(int i=0;i<str.num_nodes;i++){
		str.youngs_modulus_nodes[i] = str.youngs_modulus;
	}
}

void exePostProcess(Sim& sim, Str& str, AdjMatrix& adj_mat){
	// calculate strain and stress
	calcStressStrain(sim, str, adj_mat);
	outputStrainVtkFile(1, sim, str);
	outputStressVtkFile(1, sim, str);
	updatePosition(str);
	outputDispVtkFile(1,sim, str);
}

void initStructureStatus(Sim& sim, Str& str, AdjMatrix& adj_mat){
	// set initial state and memory
    initField(sim, str);
    initAdjMatrix(sim, str, adj_mat);
}

void exeStaticAnalysis(Sim& sim, Str& str, AdjMatrix& adj_mat){
	// output initial state
	outputParameterDataFile(sim);
	outputDispVtkFile(0, sim, str);
	outputStrainVtkFile(0, sim, str);
	outputStressVtkFile(0, sim, str);
	// calculate boundary shape function
	calcBoundaryShapeFunction(sim, str);
	// coloring Elements for Element Matrix Parallelization
	coloringElements(sim, str);
	
#ifdef MEASURE
	std::ofstream ofs(sim.time_output_filename, std::ios::app);
#endif
	// initialize sparse matrixs
	initSparseMatrix(sim, str);
	if(sim.dim == 2){
		calcElementMatrix2Dquad(sim, str, adj_mat);
	}else if(sim.dim == 3){
		calcElementMatrix3Dtetra(sim, str, adj_mat);
	}else{
		std::cout<<"Dimension Error: "<<sim.dim<<std::endl;
		exit(1);
	}
	// set sparse matrix
	setSparseMatrix(sim, str, adj_mat);
	// solve
	solveLinearEquation2D(sim, str);

#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}


void exeBucklingAnalysis(Sim& sim, Str& str, AdjMatrix& adj_mat){
	exeStaticAnalysis(sim, str, adj_mat);
	calcStressStrain(sim, str, adj_mat);
	// set stress stiffness matrix
	if(sim.dim == 2){
		calcElementMatrix2Dquad(sim, str, adj_mat);
	}else if(sim.dim == 3){
		calcElementMatrix3Dtetra(sim, str, adj_mat);
	}else{
		std::cout<<"Dimension Error: "<<sim.dim<<std::endl;
		exit(1);
	}
	// set sparse matrix
	freeSparseMatrix();
	initSparseMatrix(sim, str);
	setSparseMatrix(sim, str, adj_mat);
	// solve
	solveBuckling2D(sim, str);
	// output
	outputBucklingVtkFile(1, sim, str);
	outputStrainVtkFile(1, sim, str);
	outputStressVtkFile(1, sim, str);
	updatePosition(str);
	outputDispVtkFile(1,sim, str);
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}