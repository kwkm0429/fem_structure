#include <iostream>
#include <fstream>

#include "structure_solver.h"
#include "time_measure.h"
#include "parameter.h"

/**
 * @brief      { main }
 *
 * @return     { 0 }
 */
int main(){
#ifdef MEASURE
	initializeTime();
#endif

	solve();
	
#ifdef MEASURE
	std::ofstream ofs(sim_prm.time_output_filename, std::ios::app);
	ofs<<"SIMULATION-ENTIRE-TIME: "<<elapsedTime()/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
	return 0;
}