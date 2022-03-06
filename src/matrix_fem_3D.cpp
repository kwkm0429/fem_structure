#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>
#include <omp.h>

#include "matrix_fem_3D.h"
#include "init.h"
#include "time_measure.h"
#include "matrix_calc.h"
#include "debug.h"

void calcStiffnessMatrix3D(std::vector<std::vector<double>>& K, std::vector<std::vector<double>>& B, std::vector<std::vector<double>>& D){
	int i, j;
	std::vector<std::vector<double>> temp = std::vector< std::vector<double> >(12,std::vector<double>(3,0));
	std::vector<std::vector<double>> Bt = std::vector< std::vector<double> >(12,std::vector<double>(3,0));

	transpose_mat(Bt, B);
	multi_mat_mat(temp, Bt, D);
	multi_mat_mat(K, temp, B);
}

void calcDMatrix3DPlaneStress(std::vector<std::vector<double>>& D, Str& str){
	int i, j;
	D = {
		{1, str.poisson_ratio, 0},
		{str.poisson_ratio, 1, 0},
		{0, 0, (1-str.poisson_ratio)/2}};
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			D[i][j] *= 1 / (1 - str.poisson_ratio * str.poisson_ratio); // without Youngs modulus
		}
	}
}

void calcBMatrix3D(std::vector<std::vector<double>>& B, std::vector<double>& dN_dx, std::vector<double>& dN_dy, std::vector<double>& dN_dz, int j){
	int i;
	for(i=0;i<4;i++){
		B[0][i*3]   = dN_dx[j*4+i];
		B[1][i*3+1] = dN_dy[j*4+i];
		B[2][i*3+2] = dN_dz[j*4+i];
		B[2][i*3]   = dN_dz[j*4+i];
		B[2][i*3+1] = dN_dy[j*4+i];
		B[2][i*3+2] = dN_dx[j*4+i];
	}
}

void calcDiffMatrix3D(std::vector<std::vector<double>>& Diff, std::vector<double>& dN_dx, std::vector<double>& dN_dy, std::vector<double>& dN_dz, int j){
	int i;
	for(i=0;i<4;i++){
		Diff[0][i*3]   = dN_dx[j*4+i];
		Diff[1][i*3]   = dN_dy[j*4+i];
		Diff[2][i*3]   = dN_dz[j*4+i];
		Diff[3][i*2+1] = dN_dx[j*4+i];
		Diff[4][i*2+1] = dN_dy[j*4+i];
		Diff[5][i*2+1] = dN_dz[j*4+i];
	}
}

void calcStressMatrix3D(std::vector<std::vector<double>>& S, double sigma_x, double sigma_y, double sigma_xy){
	S[0][0] = sigma_x; S[0][1] = sigma_xy;S[0][2] = 0;       S[0][3] = 0;
	S[1][0] = sigma_xy;S[1][1] = sigma_y; S[1][2] = 0;       S[1][3] = 0;
	S[2][0] = 0;       S[2][1] = 0;       S[2][2] = sigma_x; S[2][3] = sigma_xy;
	S[3][0] = 0;       S[3][1] = 0;       S[3][2] = sigma_xy;S[3][3] = sigma_y;
}

void calcStiffGeoMatrix3D(std::vector<std::vector<double>>& Kg, std::vector<std::vector<double>>& Diff, std::vector<std::vector<double>>& S){
	int i, j;
	std::vector<std::vector<double>> temp = std::vector< std::vector<double> >(12,std::vector<double>(4,0));
	std::vector<std::vector<double>> Difft = std::vector< std::vector<double> >(12,std::vector<double>(4,0));

	transpose_mat(Difft, Diff);
	multi_mat_mat(temp, Difft, S);
	multi_mat_mat(Kg, temp, Diff);
}

void calcMassMatrix3D(std::vector<std::vector<double>>& M, std::vector<double>& N, int j, Str& str){
	int i, ii;
	for(i=0;i<4;i++){
		for(ii=0;ii<4;ii++){
			M[i*2][ii*2]     += N[j*4+i]*N[j*4+ii]*str.density;
			M[i*2][ii*2+1]   += N[j*4+i]*N[j*4+ii]*str.density;
			M[i*2+1][ii*2]   += N[j*4+i]*N[j*4+ii]*str.density;
			M[i*2+1][ii*2+1] += N[j*4+i]*N[j*4+ii]*str.density;
		}
	}
}

void calcElementMatrix3Dtetra(Sim& sim, Str& str, AdjMatrix& adj_mat){
#ifdef MEASURE
	double t_start = elapsedTime();
#endif
	int i = 0, j = 0, k = 0, l = 0;
	int size_of_column_nonzero = 0;
	int connect_check = 0;
	// tetrahedron unique parameters
	double volume = 0;
	std::vector<int> node_id = std::vector<int>(4,0);
	std::vector<std::vector<double>> mat3d = std::vector< std::vector<double> >(3,std::vector<double>(3,0));
	std::vector<std::vector<double>> mat4d = std::vector< std::vector<double> >(4,std::vector<double>(4,0));
	std::vector<double> coeff   = std::vector<double>(4, 0);
	std::vector<double> coeff_x = std::vector<double>(4, 0);
	std::vector<double> coeff_y = std::vector<double>(4, 0);
	std::vector<double> coeff_z = std::vector<double>(4, 0);
	// element matrix
	std::vector<std::vector<double>> stiff = std::vector< std::vector<double> >(4,std::vector<double>(4,0));
	// fluid local variables for OpenMP parallel
	std::vector<std::vector<double> > f_element_func = str.element_func;
	// coloring
	int elem_id = 0, color = 0;

	// initialize connectivity format
	allocateAdjMatrix(sim, str, adj_mat);

	for(color=0;color<(int)str.colored_elem_id.size();color++){
#ifdef _OPENMP
		#pragma omp parallel
		{
		#pragma omp for firstprivate( \
			j, k, l, volume, connect_check, size_of_column_nonzero, \
			node_id, mat3d, mat4d ,coeff, coeff_x, coeff_y, coeff_z, \
			stiff, \
			elem_id)
#endif 
		for(i=0;i<(int)str.colored_elem_id[color].size();i++){
			elem_id = str.colored_elem_id[color][i];
			for(j=0;j<4;j++){
				node_id[j] = str.element_node_table[elem_id][j];
			}
			// calculate volume of the element
			mat4d = {
				{1.0, 1.0, 1.0, 1.0},
				{str.x[node_id[0]], str.x[node_id[1]], str.x[node_id[2]], str.x[node_id[3]]},
				{str.y[node_id[0]], str.y[node_id[1]], str.y[node_id[2]], str.y[node_id[3]]},
				{str.z[node_id[0]], str.z[node_id[1]], str.z[node_id[2]], str.z[node_id[3]]}
			};
			volume = calc_det_4d(mat4d) / 6.0;
			if(volume<0){
				std::cout<<"Error: Negative volume"<<std::endl;
				exit(1);
			}
			// calculate stabilization parameter
			for(j=0;j<4;j++){
				f_element_func[elem_id][j] = volume / 4;
			}
			for(j=0;j<4;j++){
				// id = node_id[0]
				// lambda_id = volume_id / volume_e = a_id + b_id x_id + c_id y_id + d_id z_id
				if(j>0){
					// Change node_id order based on index of node;
					std::swap(node_id[0], node_id[1]);
					std::swap(node_id[1], node_id[2]);
					std::swap(node_id[2], node_id[3]);
				}
				// a_id
				mat3d = {
					{str.x[node_id[1]], str.x[node_id[2]], str.x[node_id[3]]},
					{str.y[node_id[1]], str.y[node_id[2]], str.y[node_id[3]]},
					{str.z[node_id[1]], str.z[node_id[2]], str.z[node_id[3]]}
				};
				if(j%2==0){
					coeff[j] =   calc_det_3d(mat3d) / 6.0 / volume;
				}else{
					coeff[j] = - calc_det_3d(mat3d) / 6.0 / volume;
				}
				// b_id
				mat3d = {
					{1.0, 1.0, 1.0},
					{str.y[node_id[1]], str.y[node_id[2]], str.y[node_id[3]]},
					{str.z[node_id[1]], str.z[node_id[2]], str.z[node_id[3]]}
				};
				if(j%2==0){
					coeff_x[j] = - calc_det_3d(mat3d) / 6.0 / volume;
				}else{
					coeff_x[j] =   calc_det_3d(mat3d) / 6.0 / volume;
				}
				// c_id
				mat3d = {
					{1.0, 1.0, 1.0},
					{str.x[node_id[1]], str.x[node_id[2]], str.x[node_id[3]]},
					{str.z[node_id[1]], str.z[node_id[2]], str.z[node_id[3]]}
				};
				if(j%2==0){
					coeff_y[j] =   calc_det_3d(mat3d) / 6.0 / volume;
				}else{
					coeff_y[j] = - calc_det_3d(mat3d) / 6.0 / volume;
				}
				// d_id
				mat3d = {
					{1.0, 1.0, 1.0},
					{str.x[node_id[1]], str.x[node_id[2]], str.x[node_id[3]]},
					{str.y[node_id[1]], str.y[node_id[2]], str.y[node_id[3]]}
				};
				if(j%2==0){
					coeff_z[j] = - calc_det_3d(mat3d) / 6.0 / volume;
				}else{
					coeff_z[j] =   calc_det_3d(mat3d) / 6.0 / volume;
				}
			}
			for(j=0;j<4;j++){
				// reset node_id
				node_id[j] = str.element_node_table[elem_id][j];
			}
			// set element matrix to adjacency matrix format
			for(j=0;j<4;j++){
				connect_check = 0;
				size_of_column_nonzero = adj_mat.idx[node_id[j]].size();
				for(k=0;k<4;k++){
					for(l=0;l<size_of_column_nonzero;l++){
						if(adj_mat.idx[node_id[j]][l] == node_id[k]){
							connect_check++;
							adj_mat.stiff[node_id[j]][l] += stiff[j][k];
							break;
						}
					}
				}
				if(connect_check != 4){
					std::cerr<<"connectivity error"<<std::endl;
					exit(1);
				}
			}
		}
#ifdef _OPENMP
		}
#endif
	}
	str.element_func = f_element_func;
#ifdef MEASURE
	double t_end = elapsedTime();
	writeTime(sim.time_output_filename, __func__, t_start, t_end);
#endif
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}