#include <iostream>
#include <fstream>

#include "structure_solver.h"
#include "time_measure.h"
#include "parameter.h"
#include "topopt.h"

int main(){
#ifdef MEASURE
	double t_start = elapsedTime();
#endif

	Sim sim; // struct for simulation setup parameter
	Str str; // struct for structure parameter
	AdjMatrix adj_mat; // adjacency matrix for connectivity
	TopOpt top; // struct for topology optimization

	readInputFiles(sim, str);
	initStructureStatus(sim, str, adj_mat);
	switch(sim.id){
		case 0: // topology optimization
			exeTopOpt(top, sim, str, adj_mat);
			break;
		case 1: // static analysis
			exeStaticAnalysis(sim, str, adj_mat);
			exePostProcess(sim, str, adj_mat);
			break;
		case 2: // buckling analysis
			exeBucklingAnalysis(sim, str, adj_mat);
			break;
		case 3: // modal analysis
			exeModalAnalysis(sim, str, adj_mat);
		default:
			break;
	}

#ifdef MEASURE
	double t_end = elapsedTime();
	writeTime(sim.time_output_filename, __func__, t_start, t_end);
#endif
	return 0;
}