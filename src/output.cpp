#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <stdlib.h> 

#include "parameter.h"
#include "output.h"
#include "matrix_calc.h"
//#include "debug.h"
#include "topopt.h"

void outputDispVtkFile(int number, Sim& sim, Str& str){
	int i,j;
	char head[16]="disp";
	char end[16]=".vtk";
	char filename[256];
	sprintf(filename,"%s%d%s",head,number,end);
	std::ofstream ofs(sim.output_data_dirname + filename);
	ofs<<"# vtk DataFile Version 4.0"<<std::endl;
	ofs<<"scalar"<<std::endl;
	ofs<<"ASCII"<<std::endl;
	ofs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	ofs<<"POINTS "<<str.num_nodes<<" float"<<std::endl;
	for(i=0;i<str.num_nodes;i++){
		ofs<<str.x[i]<<" "<<str.y[i]<<" "<<str.z[i]<<std::endl;
	}
	ofs<<"CELLS "<<str.num_elements<<" "<<5*str.num_elements<<std::endl;
	for(i=0;i<str.num_elements;i++){
		ofs<<sim.num_polygon_corner<<" ";
		for(j=0;j<sim.num_polygon_corner;j++){
			ofs<<str.element_node_table[i][j];
			if(j!=sim.num_polygon_corner-1)ofs<<" ";
			else ofs<<std::endl;
		}
	}
	ofs<<"CELL_TYPES "<<str.num_elements<<std::endl;
	for(i=0;i<str.num_elements;i++){
		if(sim.dim == 2)ofs<<9<<std::endl;
		else if(sim.dim == 3)ofs<<10<<std::endl;
	}
	ofs<<"POINT_DATA "<<str.num_nodes<<std::endl;
	ofs<<"VECTORS Displacement_Vector float"<<std::endl;
	for(i=0;i<str.num_nodes;i++){
		ofs<<str.disp_x[i]<<" "<<str.disp_y[i]<<" "<<str.disp_z[i]<<std::endl;
	}
	ofs<<"SCALARS Displacement_Scalar float"<<std::endl;
	ofs<<"LOOKUP_TABLE default"<<std::endl;
	for(i=0;i<str.num_nodes;i++){
		ofs<<std::sqrt(str.disp_x[i]*str.disp_x[i]+str.disp_y[i]*str.disp_y[i]+str.disp_z[i]*str.disp_z[i])<<std::endl;
	}
	ofs.close();
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void outputStrainVtkFile(int number, Sim& sim, Str& str){
	int i,j;
	char head[16]="strain";
	char end[16]=".vtk";
	char filename[256];
	sprintf(filename,"%s%d%s",head,number,end);
	std::ofstream ofs(sim.output_data_dirname + filename);
	ofs<<"# vtk DataFile Version 4.0"<<std::endl;
	ofs<<"scalar"<<std::endl;
	ofs<<"ASCII"<<std::endl;
	ofs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	ofs<<"POINTS "<<str.num_nodes<<" float"<<std::endl;
	for(i=0;i<str.num_nodes;i++){
		ofs<<str.x[i]<<" "<<str.y[i]<<" "<<str.z[i]<<std::endl;
	}
	ofs<<"CELLS "<<str.num_elements<<" "<<5*str.num_elements<<std::endl;
	for(i=0;i<str.num_elements;i++){
		ofs<<sim.num_polygon_corner<<" ";
		for(j=0;j<sim.num_polygon_corner;j++){
			ofs<<str.element_node_table[i][j];
			if(j!=sim.num_polygon_corner-1)ofs<<" ";
			else ofs<<std::endl;
		}
	}
	ofs<<"CELL_TYPES "<<str.num_elements<<std::endl;
	for(i=0;i<str.num_elements;i++){
		if(sim.dim == 2)ofs<<9<<std::endl;
		else if(sim.dim == 3)ofs<<10<<std::endl;
	}
	ofs<<"CELL_DATA "<<str.num_elements<<std::endl;
	ofs<<"VECTORS Strain_Vector float"<<std::endl;
	for(i=0;i<str.num_elements;i++){
		ofs<<str.strain_x[i]<<" "<<str.strain_y[i]<<" "<<str.strain_z[i]<<std::endl;
	}
	ofs<<"SCALARS Strain_Scalar float"<<std::endl;
	ofs<<"LOOKUP_TABLE default"<<std::endl;
	for(i=0;i<str.num_elements;i++){
		ofs<<std::sqrt(str.strain_x[i]*str.strain_x[i]+str.strain_y[i]*str.strain_y[i]+str.strain_z[i]*str.strain_z[i])<<std::endl;
	}
	ofs.close();
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void outputStressVtkFile(int number, Sim& sim, Str& str){
	int i,j;
	char head[16]="stress";
	char end[16]=".vtk";
	char filename[256];
	sprintf(filename,"%s%d%s",head,number,end);
	std::ofstream ofs(sim.output_data_dirname + filename);
	ofs<<"# vtk DataFile Version 4.0"<<std::endl;
	ofs<<"scalar"<<std::endl;
	ofs<<"ASCII"<<std::endl;
	ofs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	ofs<<"POINTS "<<str.num_nodes<<" float"<<std::endl;
	for(i=0;i<str.num_nodes;i++){
		ofs<<str.x[i]<<" "<<str.y[i]<<" "<<str.z[i]<<std::endl;
	}
	ofs<<"CELLS "<<str.num_elements<<" "<<5*str.num_elements<<std::endl;
	for(i=0;i<str.num_elements;i++){
		ofs<<sim.num_polygon_corner<<" ";
		for(j=0;j<sim.num_polygon_corner;j++){
			ofs<<str.element_node_table[i][j];
			if(j!=sim.num_polygon_corner-1)ofs<<" ";
			else ofs<<std::endl;
		}
	}
	ofs<<"CELL_TYPES "<<str.num_elements<<std::endl;
	for(i=0;i<str.num_elements;i++){
		if(sim.dim == 2)ofs<<9<<std::endl;
		else if(sim.dim == 3)ofs<<10<<std::endl;
	}
	ofs<<"CELL_DATA "<<str.num_elements<<std::endl;
	ofs<<"VECTORS Stress_Vector float"<<std::endl;
	for(i=0;i<str.num_elements;i++){
		ofs<<str.stress_x[i]<<" "<<str.stress_y[i]<<" "<<str.stress_z[i]<<std::endl;
	}
	ofs<<"SCALARS Stress_Scalar float"<<std::endl;
	ofs<<"LOOKUP_TABLE default"<<std::endl;
	for(i=0;i<str.num_elements;i++){
		ofs<<std::sqrt(str.stress_x[i]*str.stress_x[i]+str.stress_y[i]*str.stress_y[i]+str.stress_z[i]*str.stress_z[i])<<std::endl;
	}
	ofs.close();
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void outputBucklingVtkFile(int number, Sim& sim, Str& str){
	int i,j;
	char head[16]="buckling";
	char end[16]=".vtk";
	char filename[256];
	sprintf(filename,"%s%d%s",head,number,end);
	std::ofstream ofs(sim.output_data_dirname + filename);
	ofs<<"# vtk DataFile Version 4.0"<<std::endl;
	ofs<<"scalar"<<std::endl;
	ofs<<"ASCII"<<std::endl;
	ofs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	ofs<<"POINTS "<<str.num_nodes<<" float"<<std::endl;
	for(i=0;i<str.num_nodes;i++){
		ofs<<str.x[i]+str.buckling_x[i]<<" "<<str.y[i]+str.buckling_y[i]<<" "<<str.z[i]+str.buckling_z[i]<<std::endl;
	}
	ofs<<"CELLS "<<str.num_elements<<" "<<5*str.num_elements<<std::endl;
	for(i=0;i<str.num_elements;i++){
		ofs<<sim.num_polygon_corner<<" ";
		for(j=0;j<sim.num_polygon_corner;j++){
			ofs<<str.element_node_table[i][j];
			if(j!=sim.num_polygon_corner-1)ofs<<" ";
			else ofs<<std::endl;
		}
	}
	ofs<<"CELL_TYPES "<<str.num_elements<<std::endl;
	for(i=0;i<str.num_elements;i++){
		if(sim.dim == 2)ofs<<9<<std::endl;
		else if(sim.dim == 3)ofs<<10<<std::endl;
	}
	ofs<<"POINT_DATA "<<str.num_nodes<<std::endl;
	ofs<<"VECTORS Buckling_Vector float"<<std::endl;
	for(i=0;i<str.num_nodes;i++){
		ofs<<str.buckling_x[i]<<" "<<str.buckling_y[i]<<" "<<str.buckling_z[i]<<std::endl;
	}
	ofs<<"SCALARS Buckling_Scalar float"<<std::endl;
	ofs<<"LOOKUP_TABLE default"<<std::endl;
	for(i=0;i<str.num_nodes;i++){
		ofs<<std::sqrt(str.buckling_x[i]*str.buckling_x[i]+str.buckling_y[i]*str.buckling_y[i]+str.buckling_z[i]*str.buckling_z[i])<<std::endl;
	}
	ofs.close();
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}


void outputDensityVtkFile(int number, TopOpt& top, Sim& sim, Str& str){
	int i,j;
	char head[16]="density";
	char end[16]=".vtk";
	char filename[256];
	sprintf(filename,"%s%d%s",head,number,end);
	std::ofstream ofs(sim.output_data_dirname + filename);
	ofs<<"# vtk DataFile Version 4.0"<<std::endl;
	ofs<<"scalar"<<std::endl;
	ofs<<"ASCII"<<std::endl;
	ofs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	ofs<<"POINTS "<<str.num_nodes<<" float"<<std::endl;
	for(i=0;i<str.num_nodes;i++){
		ofs<<str.x[i]<<" "<<str.y[i]<<" "<<str.z[i]<<std::endl;
	}
	ofs<<"CELLS "<<str.num_elements<<" "<<5*str.num_elements<<std::endl;
	for(i=0;i<str.num_elements;i++){
		ofs<<sim.num_polygon_corner<<" ";
		for(j=0;j<sim.num_polygon_corner;j++){
			ofs<<str.element_node_table[i][j];
			if(j!=sim.num_polygon_corner-1)ofs<<" ";
			else ofs<<std::endl;
		}
	}
	ofs<<"CELL_TYPES "<<str.num_elements<<std::endl;
	for(i=0;i<str.num_elements;i++){
		if(sim.dim == 2)ofs<<9<<std::endl;
		else if(sim.dim == 3)ofs<<10<<std::endl;
	}
	ofs<<"POINT_DATA "<<str.num_nodes<<std::endl;
	ofs<<"SCALARS Density float"<<std::endl;
	ofs<<"LOOKUP_TABLE default"<<std::endl;
	for(i=0;i<str.num_nodes;i++){
		ofs<<top.rho[i]<<std::endl;
	}
	ofs.close();
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void outputParameterDataFile(Sim& sim){
	char filename[32]="param.log";
	std::ofstream ofs(sim.output_data_dirname + filename);
	ofs<<"DIMENSION = "<<sim.dim<<std::endl;
	ofs<<"SHAPE_NUMBER = "<<sim.num_polygon_corner<<std::endl;
	ofs<<"TIME_STEP = "<<sim.max_time_step<<std::endl;
	ofs<<"OUTPUT_INTERVAL = "<<sim.output_interval<<std::endl;
	ofs<<"DT = "<<sim.dt<<std::endl;
	ofs<<"GRAVITY_X = "<<sim.gravity_x<<std::endl;
	ofs<<"GRAVITY_Y = "<<sim.gravity_y<<std::endl;
	ofs<<"GRAVITY_Z = "<<sim.gravity_z<<std::endl;
	ofs<<"EQUATION_SOLVER = "<<sim.eq_solver_opt<<std::endl;
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

void outputTopOptDataFile(TopOpt& top, Sim& sim){
	char filename[32]="top.log";
	std::ofstream ofs(sim.output_data_dirname + filename);
	ofs<<"VOL_MAX = "<<top.vol_max<<std::endl;
	ofs<<"RHO_INIT = "<<top.rho_init<<std::endl;
	ofs<<"E0 = "<<top.E0<<std::endl;
	ofs<<"Emin = "<<top.Emin<<std::endl;
	ofs<<"POW = "<<top.pow<<std::endl;
	ofs<<"ITR_MAX = "<<top.itr_max<<std::endl;
	ofs<<"FILTER_RADIUS = "<<top.filter_radius<<std::endl;
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