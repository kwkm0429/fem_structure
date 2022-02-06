#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <stdlib.h> 

#include "parameter.h"
#include "output.h"
#include "matrix_calc.h"
#include "debug.h"

void outputDispVtkFile(int number){
	int i,j;
	char head[16]="disp";
	char end[16]=".vtk";
	char filename[256];
	sprintf(filename,"%s%d%s",head,number,end);
	std::ofstream ofs(sim_prm.output_data_dirname + filename);
	ofs<<"# vtk DataFile Version 4.0"<<std::endl;
	ofs<<"scalar"<<std::endl;
	ofs<<"ASCII"<<std::endl;
	ofs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	ofs<<"POINTS "<<structure.num_nodes<<" float"<<std::endl;
	for(i=0;i<structure.num_nodes;i++){
		ofs<<structure.x[i]<<" "<<structure.y[i]<<" "<<structure.z[i]<<std::endl;
	}
	ofs<<"CELLS "<<structure.num_elements<<" "<<5*structure.num_elements<<std::endl;
	for(i=0;i<structure.num_elements;i++){
		ofs<<sim_prm.num_polygon_corner<<" ";
		for(j=0;j<sim_prm.num_polygon_corner;j++){
			ofs<<structure.element_node_table[i][j];
			if(j!=sim_prm.num_polygon_corner-1)ofs<<" ";
			else ofs<<std::endl;
		}
	}
	ofs<<"CELL_TYPES "<<structure.num_elements<<std::endl;
	for(i=0;i<structure.num_elements;i++){
		if(sim_prm.dim == 2)ofs<<9<<std::endl;
		else if(sim_prm.dim == 3)ofs<<10<<std::endl;
	}
	ofs<<"POINT_DATA "<<structure.num_nodes<<std::endl;
	ofs<<"VECTORS velocity float"<<std::endl;
	for(i=0;i<structure.num_nodes;i++){
		ofs<<structure.disp_x[i]<<" "<<structure.disp_y[i]<<" "<<structure.disp_z[i]<<std::endl;
	}
	ofs<<"SCALARS point_scalars float"<<std::endl;
	ofs<<"LOOKUP_TABLE default"<<std::endl;
	for(i=0;i<structure.num_nodes;i++){
		ofs<<std::sqrt(structure.disp_x[i]*structure.disp_x[i]+structure.disp_y[i]*structure.disp_y[i]+structure.disp_z[i]*structure.disp_z[i])<<std::endl;
	}
	ofs.close();
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void outputParameterDataFile(){
	char filename[32]="param.log";
	std::ofstream ofs(sim_prm.output_data_dirname + filename);
	ofs<<"DIMENSION = "<<sim_prm.dim<<std::endl;
	ofs<<"SHAPE_NUMBER = "<<sim_prm.num_polygon_corner<<std::endl;
	ofs<<"TIME_STEP = "<<sim_prm.max_time_step<<std::endl;
	ofs<<"OUTPUT_INTERVAL = "<<sim_prm.output_interval<<std::endl;
	ofs<<"DT = "<<sim_prm.dt<<std::endl;
	ofs<<"GRAVITY_X = "<<sim_prm.gravity_x<<std::endl;
	ofs<<"GRAVITY_Y = "<<sim_prm.gravity_y<<std::endl;
	ofs<<"GRAVITY_Z = "<<sim_prm.gravity_z<<std::endl;
	ofs<<"FLUID_DENSITY = "<<structure.density<<std::endl;
	ofs<<"EQUATION_SOLVER = "<<sim_prm.eq_solver_opt<<std::endl;
	char* num_threads = getenv("OMP_NUM_THREADS");
	if(num_threads!=NULL){
		ofs<<"OMP_NUM_THREADS = "<<num_threads<<std::endl;
	}else{
		ofs<<"OMP_NUM_THREADS = NULL"<<std::endl;
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
	ofs.close();
}