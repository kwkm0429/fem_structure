#include <iostream>
#include <fstream>

#include "structure_solver.h"
#include "time_measure.h"
#include "parameter.h"
#include "topopt.h"

int main(){
#ifdef MEASURE
	initializeTime();
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
		default:
			break;
	}

#ifdef MEASURE
	std::ofstream ofs(sim.time_output_filename, std::ios::app);
	ofs<<"SIMULATION-ENTIRE-TIME: "<<elapsedTime()/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
	return 0;
}