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

	readInputFiles();
	if(sim_prm.is_topopt){
		exeTopOpt();
	}else{
		exeStaticAnalysis();
		exePostProcess();
	}

#ifdef MEASURE
	std::ofstream ofs(sim_prm.time_output_filename, std::ios::app);
	ofs<<"SIMULATION-ENTIRE-TIME: "<<elapsedTime()/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
	return 0;
}