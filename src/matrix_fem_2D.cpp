#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>
#include <omp.h>

#include "matrix_fem_2D.h"
#include "init.h"
#include "time_measure.h"
#include "eigen_solver.h"
#include "matrix_calc.h"
#include "debug.h"

void calcJacobian(
	int k, 
	std::vector<double>& J,
	std::vector<double>& x,
	std::vector<double>& y,
	std::vector<double>& N,
	std::vector<double>& dN_dx,
	std::vector<double>& dN_dy)
	{
	int l;
	double eta, xi, dx_dxi, dx_deta, dy_dxi, dy_deta;
	std::vector<double> dN_deta = std::vector<double>(4,0);
	std::vector<double> dN_dxi  = std::vector<double>(4,0);
	
	if(k==0 || k==3){
		xi = 1.0/std::sqrt(3);
	}else{
		xi = -1.0/std::sqrt(3);
	}
	if(k==0 || k==1){
		eta = 1.0/std::sqrt(3);
	}else{
		eta = -1.0/std::sqrt(3);
	}

	N[k*4+0] = (1 - xi) * (1 - eta) * 0.25;
	N[k*4+1] = (1 + xi) * (1 - eta) * 0.25;
	N[k*4+2] = (1 + xi) * (1 + eta) * 0.25;
	N[k*4+3] = (1 - xi) * (1 + eta) * 0.25;

	dN_dxi[0]  = -(1 - eta) * 0.25;
	dN_dxi[1]  =  (1 - eta) * 0.25;
	dN_dxi[2]  =  (1 + eta) * 0.25;
	dN_dxi[3]  = -(1 + eta) * 0.25;

	dN_deta[0] = -(1 - xi)  * 0.25;
	dN_deta[1] = -(1 + xi)  * 0.25;
	dN_deta[2] =  (1 + xi)  * 0.25;
	dN_deta[3] =  (1 - xi)  * 0.25;

	dx_dxi = dx_deta = dy_dxi = dy_deta = 0;

	for(l=0;l<4;l++){
		dx_dxi  += dN_dxi[l]  * x[l];
		dx_deta += dN_deta[l] * x[l];
		dy_dxi  += dN_dxi[l]  * y[l];
		dy_deta += dN_deta[l] * y[l];
	}

	J[k] = dx_dxi * dy_deta - dx_deta * dy_dxi; //Jacobian

	for(l=0;l<4;l++){
		dN_dx[k*4+l] = (dy_deta * dN_dxi[l] - dy_dxi * dN_deta[l])  / J[k];
		dN_dy[k*4+l] = (-dx_deta * dN_dxi[l] + dx_dxi * dN_deta[l]) / J[k];
	}
}

void calcStiffnessMatrix2D(std::vector<std::vector<double>>& K, std::vector<std::vector<double>>& B, std::vector<std::vector<double>>& D){
	int i, j;
	std::vector<std::vector<double>> temp = std::vector< std::vector<double> >(8,std::vector<double>(3,0));
	std::vector<std::vector<double>> Bt = std::vector< std::vector<double> >(8,std::vector<double>(3,0));

	transpose_mat(Bt, B);
	multi_mat_mat(temp, Bt, D);
	multi_mat_mat(K, temp, B);
}

void calcDMatrix2DPlaneStress(std::vector<std::vector<double>>& D, Str& str){
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

void calcBMatrix2D(std::vector<std::vector<double>>& B, std::vector<double>& dN_dx, std::vector<double>& dN_dy, int j){
	int i;
	for(i=0;i<4;i++){
		B[0][i*2]   = dN_dx[j*4+i];
		B[1][i*2+1] = dN_dy[j*4+i];
		B[2][i*2]   = dN_dy[j*4+i];
		B[2][i*2+1] = dN_dx[j*4+i];
	}
}

void calcDiffMatrix2D(std::vector<std::vector<double>>& Diff, std::vector<double>& dN_dx, std::vector<double>& dN_dy, int j){
	int i;
	for(i=0;i<4;i++){
		Diff[0][i*2]   = dN_dx[j*4+i];
		Diff[1][i*2]   = dN_dy[j*4+i];
		Diff[2][i*2+1] = dN_dx[j*4+i];
		Diff[3][i*2+1] = dN_dy[j*4+i];
	}
}

void calcStressMatrix2D(std::vector<std::vector<double>>& S, double sigma_x, double sigma_y, double sigma_xy){
	S[0][0] = sigma_x; S[0][1] = sigma_xy;S[0][2] = 0;       S[0][3] = 0;
	S[1][0] = sigma_xy;S[1][1] = sigma_y; S[1][2] = 0;       S[1][3] = 0;
	S[2][0] = 0;       S[2][1] = 0;       S[2][2] = sigma_x; S[2][3] = sigma_xy;
	S[3][0] = 0;       S[3][1] = 0;       S[3][2] = sigma_xy;S[3][3] = sigma_y;
}

void calcStiffGeoMatrix2D(std::vector<std::vector<double>>& Kg, std::vector<std::vector<double>>& Diff, std::vector<std::vector<double>>& S){
	int i, j;
	std::vector<std::vector<double>> temp = std::vector< std::vector<double> >(8,std::vector<double>(4,0));
	std::vector<std::vector<double>> Difft = std::vector< std::vector<double> >(8,std::vector<double>(4,0));

	transpose_mat(Difft, Diff);
	multi_mat_mat(temp, Difft, S);
	multi_mat_mat(Kg, temp, Diff);
}

void calcMassMatrix2D(std::vector<std::vector<double>>& M, std::vector<double>& N, int j, Str& str){
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

void calcStressStrain(Sim& sim, Str& str, AdjMatrix& adj_mat){
#ifdef MEASURE
	double t_start = elapsedTime();
#endif
	int i = 0, j = 0;
	std::vector<int> node_id = std::vector<int>(4,0);
	std::vector<double> x     = std::vector<double>(4,0);
	std::vector<double> y     = std::vector<double>(4,0);
	std::vector<double> dN_dx = std::vector<double>(16,0);
	std::vector<double> dN_dy = std::vector<double>(16,0);
	std::vector<double> N     = std::vector<double>(16,0);
	std::vector<double> J     = std::vector<double>(4,0);
	std::vector<double> disp_elem   = std::vector<double>(8,0);
	std::vector<double> strain_elem = std::vector<double>(3,0);
	std::vector<double> stress_elem = std::vector<double>(3,0);
	std::vector<std::vector<double>> strain_disp_matrix = std::vector< std::vector<double> >(3,std::vector<double>(8,0)); // B matrix
	std::vector<std::vector<double>> stress_strain_matrix = std::vector< std::vector<double> >(3,std::vector<double>(3,0));	// D matrix
	// coloring
	int elem_id = 0, color = 0;

    // calculate element stiffness matrix
	for(color=0;color<(int)str.colored_elem_id.size();color++){
#ifdef _OPENMP
		#pragma omp parallel
		{
		#pragma omp for firstprivate(node_id, j, x, y, dN_dx, dN_dy, N, J, \
		elem_id, strain_disp_matrix, stress_strain_matrix, disp_elem, strain_elem, stress_elem)
#endif 
		for(i=0;i<(int)str.colored_elem_id[color].size();i++){
			elem_id = str.colored_elem_id[color][i];
			// calculate stabilization parameter in element i
			for(j=0;j<4;j++){
				node_id[j] = str.element_node_table[elem_id][j];
				x[j] = str.x[node_id[j]];
				y[j] = str.y[node_id[j]];
				disp_elem[j*2]   = str.disp_x[node_id[j]];
				disp_elem[j*2+1] = str.disp_y[node_id[j]];
			}
			for(j=0;j<4;j++){
				calcJacobian(j,J,x,y,N,dN_dx,dN_dy); // calculate jaccobian and basis function
			}
			// calc D matrix (plane stress)
			calcDMatrix2DPlaneStress(stress_strain_matrix, str);
			// calc B matrix
			for(j=0;j<4;j++){ // quadrature point
				calcBMatrix2D(strain_disp_matrix, dN_dx, dN_dy, j);
				// calc strain and stress
				multi_mat_vec(strain_elem, strain_disp_matrix, disp_elem);
				multi_mat_vec(stress_elem, stress_strain_matrix, strain_elem);
				str.strain_x[elem_id] += strain_elem[0] * J[j] * str.youngs_modulus_nodes[node_id[j]];
				str.strain_y[elem_id] += strain_elem[1] * J[j] * str.youngs_modulus_nodes[node_id[j]];
				str.sheer_strain_xy[elem_id] += strain_elem[2] * J[j] * str.youngs_modulus_nodes[node_id[j]];
				str.stress_x[elem_id] += stress_elem[0] * J[j] * str.youngs_modulus_nodes[node_id[j]];
				str.stress_y[elem_id] += stress_elem[1] * J[j] * str.youngs_modulus_nodes[node_id[j]];
				str.sheer_stress_xy[elem_id] += stress_elem[2] * J[j] * str.youngs_modulus_nodes[node_id[j]];
			}
		}
#ifdef _OPENMP
		}
#endif
	}
#ifdef MEASURE
	double t_end = elapsedTime();
	writeTime(sim.time_output_filename, __func__, t_start, t_end);
#endif
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void calcElementMatrix2Dquad(Sim& sim, Str& str, AdjMatrix& adj_mat){
#ifdef MEASURE
	double t_start = elapsedTime();
#endif
	int i = 0, j = 0, k = 0, l = 0;
	std::vector<int> node_id = std::vector<int>(4,0);
	std::vector<double> x     = std::vector<double>(4,0);
	std::vector<double> y     = std::vector<double>(4,0);
	std::vector<double> dN_dx = std::vector<double>(16,0);
	std::vector<double> dN_dy = std::vector<double>(16,0);
	std::vector<double> N     = std::vector<double>(16,0);
	std::vector<double> J     = std::vector<double>(4,0);
	std::vector<double> disp_elem   = std::vector<double>(8,0);
	std::vector<double> strain_elem = std::vector<double>(3,0);
	std::vector<double> stress_elem = std::vector<double>(3,0);
	std::vector<std::vector<double>> strain_disp_matrix = std::vector< std::vector<double> >(3,std::vector<double>(8,0)); // B matrix
	std::vector<std::vector<double>> stress_strain_matrix = std::vector< std::vector<double> >(3,std::vector<double>(3,0));	// D matrix
	std::vector<std::vector<double>> stiff_matrix = std::vector< std::vector<double> >(8,std::vector<double>(8,0)); // K matrix
	std::vector<std::vector<double>> stiff_matrix_temp = std::vector< std::vector<double> >(8,std::vector<double>(8,0)); // K temp
	// topology optimization
	std::vector<std::vector<double>> stiff_top_matrix = std::vector< std::vector<double> >(8,std::vector<double>(8,0));
	// Kg matrix
	std::vector<std::vector<double>> stiff_geo_matrix = std::vector< std::vector<double> >(8,std::vector<double>(8,0)); // Kg matrix
	std::vector<std::vector<double>> stiff_geo_matrix_temp = std::vector< std::vector<double> >(8,std::vector<double>(8,0)); // Kg matrix
	std::vector<std::vector<double>> stress_matrix = std::vector< std::vector<double> >(4,std::vector<double>(4,0)); // stress matrix
	std::vector<std::vector<double>> diff_matrix = std::vector< std::vector<double> >(4,std::vector<double>(8,0)); // differential matrix
	// M matrix
	std::vector<std::vector<double>> mass_matrix = std::vector< std::vector<double> >(8,std::vector<double>(8,0)); // M matrix
	int connect_check = 0;
	int size_of_column_nonzero = 0;
	// fluid local variables for OpenMP parallel
	std::vector<std::vector<double> > f_element_func = str.element_func;
	// coloring
	int elem_id = 0, color = 0;
	// topology optimization
	std::vector<double> temp = std::vector<double>(8,0);
	std::vector<double> disp_top = std::vector<double>(8,0);
	double sens;

	// initialize connectivity format
	allocateAdjMatrix(sim, str, adj_mat);

    // calculate element stiffness matrix
	for(color=0;color<(int)str.colored_elem_id.size();color++){
#ifdef _OPENMP
		#pragma omp parallel
		{
		#pragma omp for firstprivate(node_id, j, k, l, \
		x, y, dN_dx, dN_dy, N, J, detJ, \
		connect_check, size_of_column_nonzero,\
		elem_id, strain_disp_matrix, stress_strain_matrix, stiff_matrix, disp_elem, strain_elem, stress_elem, temp, sens)
#endif 
		for(i=0;i<(int)str.colored_elem_id[color].size();i++){
			elem_id = str.colored_elem_id[color][i];
			// calculate stabilization parameter in element i
			for(j=0;j<4;j++){
				node_id[j] = str.element_node_table[elem_id][j];
				x[j] = str.x[node_id[j]];
				y[j] = str.y[node_id[j]];
				disp_elem[j*2]   = str.disp_x[node_id[j]];
				disp_elem[j*2+1] = str.disp_y[node_id[j]];
			}
			for(j=0;j<4;j++){
				calcJacobian(j,J,x,y,N,dN_dx,dN_dy); // calculate jaccobian and basis function
			}
			for(j=0;j<4;j++){ // nodes in the element
				f_element_func[elem_id][j] = 0;
				for(k=0;k<4;k++){ // quadrature point
					// N[l*4+k] : weight at quadrature point l in terms of node k
					f_element_func[elem_id][j] += N[k*4+j] * std::abs(J[k]);
				}
			}
			stiff_matrix = std::vector< std::vector<double> >(8,std::vector<double>(8,0)); // K matrix
			stiff_geo_matrix = std::vector< std::vector<double> >(8,std::vector<double>(8,0)); // Kg matrix
			mass_matrix = std::vector< std::vector<double> >(8,std::vector<double>(8,0)); // M matrix
			// calc D matrix (plane stress)
			calcDMatrix2DPlaneStress(stress_strain_matrix, str);
			// calc B matrix
			for(j=0;j<4;j++){ // quadrature point
				// calc M matrix
				calcMassMatrix2D(mass_matrix, N, j, str);
				// calc B matrix
				calcBMatrix2D(strain_disp_matrix, dN_dx, dN_dy, j);
				// calculate element stiffness matrix
				calcStiffnessMatrix2D(stiff_matrix_temp, strain_disp_matrix, stress_strain_matrix);
				for(k=0;k<stiff_matrix.size();k++){
					for(l=0;l<stiff_matrix[k].size();l++){
						stiff_matrix[k][l] += stiff_matrix_temp[k][l] * str.thickness * J[j] * str.youngs_modulus_nodes[node_id[j]];
						stiff_top_matrix[k][l] = stiff_matrix_temp[k][l] * str.thickness * J[j];
					}
				}
				// calc stiff_geo matrix
				calcDiffMatrix2D(diff_matrix, dN_dx, dN_dy, j);
				calcStressMatrix2D(stress_matrix, str.stress_x[elem_id], str.stress_y[elem_id], str.sheer_stress_xy[elem_id]);
				calcStiffGeoMatrix2D(stiff_geo_matrix_temp, diff_matrix, stress_matrix);
				for(k=0;k<stiff_geo_matrix.size();k++){
					for(l=0;l<stiff_geo_matrix.size();l++){
						stiff_geo_matrix[k][l] += stiff_geo_matrix_temp[k][l] * str.thickness * J[j];
					}
				}
				// calc topology optimization sensitivity
				multi_mat_vec(temp, stiff_top_matrix, disp_elem);
				sens = multi_vec_vec(disp_elem, temp);
				for(k=0;k<4;k++){
					str.sensitivity[node_id[k]] = sens;
				}
			}
			// set element matrix to adjacency matrix format
			for(j=0;j<4;j++){
				connect_check = 0;
				size_of_column_nonzero = adj_mat.idx[node_id[j]].size();
				for(k=0;k<4;k++){
					for(l=0;l<size_of_column_nonzero;l++){
						if(adj_mat.idx[node_id[j]][l] == node_id[k]){
							connect_check++;
							// linear stiffness matrix
							adj_mat.stiff[node_id[j]*2][l*2]     += stiff_matrix[j*2][k*2];
							adj_mat.stiff[node_id[j]*2][l*2+1]   += stiff_matrix[j*2][k*2+1];
							adj_mat.stiff[node_id[j]*2+1][l*2]   += stiff_matrix[j*2+1][k*2];
							adj_mat.stiff[node_id[j]*2+1][l*2+1] += stiff_matrix[j*2+1][k*2+1];
							// geometry stiffness matrix
							adj_mat.stiff_geo[node_id[j]*2][l*2]     += stiff_geo_matrix[j*2][k*2];
							adj_mat.stiff_geo[node_id[j]*2][l*2+1]   += stiff_geo_matrix[j*2][k*2+1];
							adj_mat.stiff_geo[node_id[j]*2+1][l*2]   += stiff_geo_matrix[j*2+1][k*2];
							adj_mat.stiff_geo[node_id[j]*2+1][l*2+1] += stiff_geo_matrix[j*2+1][k*2+1];
							// mass matrix
							adj_mat.mass[node_id[j]*2][l*2]     += mass_matrix[j*2][k*2];
							adj_mat.mass[node_id[j]*2][l*2+1]   += mass_matrix[j*2][k*2+1];
							adj_mat.mass[node_id[j]*2+1][l*2]   += mass_matrix[j*2+1][k*2];
							adj_mat.mass[node_id[j]*2+1][l*2+1] += mass_matrix[j*2+1][k*2+1];
						}
					}
				}
				if(connect_check != 4){
					std::cerr<<"connectivity error: "<<connect_check<<std::endl;
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
	std::cout<<"Matrix nonzero element num : "<<sim.num_nonzero<<std::endl;
	std::cout<<"Matrix   total element num : "<<(long long)(str.num_nodes)*str.num_nodes<<"\n"<<std::endl;
    debugPrintInfo(__func__);
#endif
}

void calcBoundaryShapeFunction(Sim& sim, Str& str){
	int i = 0, j = 0, node_id;
	std::vector<double> pos_x, pos_y, pos_z;
	double dx1, dy1, dz1, dx2, dy2, dz2, dd1, dd2, d1d2, area = 0;

	// fluid local variables for OpenMP parallel
	for(i=0;i<str.num_elements;i++){
		// check if element i include boundary surface
		for(j=0;j<4;j++){
			node_id = str.element_node_table[i][j];
			if(str.is_boundary[node_id]){
				pos_x.push_back(str.x[node_id]);
				pos_y.push_back(str.y[node_id]);
				pos_z.push_back(str.z[node_id]);
			}
		}
		if(sim.dim == 2 && pos_x.size() == 2){
			dx1 = pos_x[1] - pos_x[0];
			dy1 = pos_y[1] - pos_y[0];
			dz1 = pos_z[1] - pos_z[0];
			area = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1); // area means length in 2D
			for(j=0;j<4;j++){
				node_id = str.element_node_table[i][j];
				if(str.is_boundary[node_id]){
					str.boundary_shape_function[node_id] += area / 2.0;
				}
			}
		}else if(sim.dim == 3 && pos_x.size() == 3){
			dx1 = pos_x[1] - pos_x[0];
			dy1 = pos_y[1] - pos_y[0];
			dz1 = pos_z[1] - pos_z[0];
			dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
			dx2 = pos_x[2] - pos_x[0];
			dy2 = pos_y[2] - pos_y[0];
			dz2 = pos_z[2] - pos_z[0];
			dd2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
			d1d2 = dx1*dx2 + dy1*dy2 + dz1*dz2;
			area = 0.5 * std::sqrt(dd1*dd2 - d1d2*d1d2); // area of a triangle in 3D
			for(j=0;j<4;j++){
				node_id = str.element_node_table[i][j];
				if(str.is_boundary[node_id]){
					str.boundary_shape_function[node_id] += area / 3.0;
				}
			}
		}
		pos_x.clear();
		pos_y.clear();
		pos_z.clear();
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

/**
 * @brief   Coloring elements for parallelization. Elements are colored by Welsh-Powell Algorithm.
 */
void coloringElements(Sim& sim, Str& str){
	std::vector<std::vector<int>> node2elem = std::vector<std::vector<int>>(str.num_nodes);
	std::vector<std::vector<int>> elem2node = std::vector<std::vector<int>>(str.num_elements);
	std::vector<std::vector<int>> elem2elem = std::vector<std::vector<int>>(str.num_elements);
	std::vector<std::pair<int, int>> deg_size = std::vector<std::pair<int, int>>(str.num_elements);
	std::vector<int> color_elem = std::vector<int>(str.num_elements, -1);
	int i, j, k, node_id, elem_id1, elem_id2, max_color = 0;

	// create node-element table
	for(i=0;i<str.num_elements;i++){
		for(j=0;j<sim.num_polygon_corner;j++){
			node_id = str.element_node_table[i][j];
			node2elem[node_id].push_back(i);
		}
	}
	// create elemet-node table
	for(i=0;i<str.num_nodes;i++){
		for(j=0;j<(int)node2elem[i].size();j++){
			elem_id1 = node2elem[i][j];
			elem2node[elem_id1].push_back(i);
		}
	}
	// create element-element table (graph)
	for(i=0;i<str.num_elements;i++){
		std::vector<bool> is_registered = std::vector<bool>(str.num_elements, false);
		for(j=0;j<(int)elem2node[i].size();j++){
			node_id = elem2node[i][j];
			for(k=0;k<(int)node2elem[node_id].size();k++){
				elem_id1 = node2elem[node_id][k];
				if(i==elem_id1)continue;
				if(is_registered[elem_id1])continue;
				elem2elem[i].push_back(elem_id1);
				is_registered[elem_id1] = true;
			}
		}
	}
	// sort by number of nodes belong to the element
	for(i=0;i<str.num_elements;i++){
		deg_size[i] = std::pair<int, int>(elem2elem[i].size(), i);
	}
	sort(deg_size.rbegin(), deg_size.rend());
	// Welsh-Powell Algorithm
	for(i=0;i<(int)deg_size.size();i++){
		elem_id1 = deg_size[i].second;
		int color = 0;
		bool is_used = true;
		while(is_used){
			is_used = false;
			for(j=0;j<(int)elem2elem[elem_id1].size();j++){
				elem_id2 = elem2elem[elem_id1][j];
				if(elem_id1 == elem_id2)continue;
				if(color_elem[elem_id2] == color){
					is_used = true;
					color++;
					break;
				}
			}
		}
		color_elem[elem_id1] = color;
		max_color = std::max(max_color, color);
	}
	// set element group by color
	str.colored_elem_id = std::vector<std::vector<int>>(max_color+1);
	for(i=0;i<str.num_elements;i++){
		int color = color_elem[i];
		str.colored_elem_id[color].push_back(i);
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
    std::cout<<"Number of colors: "<<max_color<<std::endl;
	for(i=0;i<max_color;i++){
		std::cout<<"Color id: "<<i<<" "<<"Number of nodes: "<<str.colored_elem_id[i].size()<<std::endl;
	}
#endif
}