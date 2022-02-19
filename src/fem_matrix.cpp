#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>
#include <omp.h>

#include "parameter.h"
#include "fem_matrix.h"
#include "init.h"
#include "time_measure.h"
#include "eigen_solver.h"
#include "matrix_calc.h"
//#include "debug.h"

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

void calcStiffnessMatrix(
	std::vector<std::vector<double>>& K,
	std::vector<std::vector<double>>& B,
	std::vector<std::vector<double>>& D){
	int i, j;
	std::vector<std::vector<double>> temp = std::vector< std::vector<double> >(8,std::vector<double>(3,0));
	std::vector<std::vector<double>> Bt = std::vector< std::vector<double> >(8,std::vector<double>(3,0));
	for(i=0;i<8;i++){
		for(j=0;j<3;j++){
			Bt[i][j] = B[j][i];
		}
	}
	multi_mat_mat(temp, Bt, D);
	multi_mat_mat(K, temp, B);
}

void calcElementMatrix2Dquad(){
#ifdef MEASURE
	double time_start = elapsedTime();
#endif
	int node_id1 = 0, node_id2 = 0, node_id3 = 0, node_id4 = 0;
	int i = 0, j = 0, k = 0, l = 0, m=0, ne1 = 0, ne2 = 0, ne3 = 0, ne4 = 0;
	double detJ = 0;
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
	std::vector<std::vector<double>> stiff_top_matrix = std::vector< std::vector<double> >(8,std::vector<double>(8,0));

	double mass = 0, stiff = 0, damping = 0;
	int connect_check = 0;
	int size_of_column_nonzero = 0;
	// fluid local variables for OpenMP parallel
	std::vector<std::vector<double> > f_element_func = structure.element_func;
	// coloring
	int elem_id = 0, color = 0;
	// topology optimization
	std::vector<double> temp = std::vector<double>(8,0);
	std::vector<double> disp_top = std::vector<double>(8,0);
	double sens;

	// initialize connectivity format
	allocateAdjMatrix();

    // calculate element stiffness matrix
	for(color=0;color<(int)structure.colored_elem_id.size();color++){
#ifdef _OPENMP
		#pragma omp parallel
		{
		#pragma omp for firstprivate(node_id, node_id1, node_id2, node_id3, node_id4, \
		ne1, ne2, ne3, ne4, j, k, l, m, \
		x, y, dN_dx, dN_dy, N, J, detJ, \
		mass, stiff, damping, connect_check, size_of_column_nonzero,\
		elem_id, strain_disp_matrix, stress_strain_matrix, stiff_matrix, disp_elem, strain_elem, stress_elem, temp, sens)
#endif 
		for(i=0;i<(int)structure.colored_elem_id[color].size();i++){
			elem_id = structure.colored_elem_id[color][i];
			// calculate stabilization parameter in element i
			for(j=0;j<4;j++){
				node_id[j] = structure.element_node_table[elem_id][j];
				x[j] = structure.x[node_id[j]];
				y[j] = structure.y[node_id[j]];
				disp_elem[j*2]   = structure.disp_x[node_id[j]];
				disp_elem[j*2+1] = structure.disp_y[node_id[j]];
			}
			detJ=0;
			for(j=0;j<4;j++){
				calcJacobian(j,J,x,y,N,dN_dx,dN_dy); // calculate jaccobian and basis function
				detJ += std::abs(J[j]);
			}
			for(j=0;j<4;j++){ // nodes in the element
				f_element_func[elem_id][j] = 0;
				for(k=0;k<4;k++){ // quadrature point
					// N[l*4+k] : weight at quadrature point l in terms of node k
					f_element_func[elem_id][j] += N[k*4+j] * std::abs(J[k]);
				}
			}
			for(j=0;j<4;j++){ // quadrature point
				strain_disp_matrix[0][0] = dN_dx[j*4+0];
				strain_disp_matrix[0][2] = dN_dx[j*4+1];
				strain_disp_matrix[0][4] = dN_dx[j*4+2];
				strain_disp_matrix[0][6] = dN_dx[j*4+3];
				
				strain_disp_matrix[1][1] = dN_dy[j*4+0];
				strain_disp_matrix[1][3] = dN_dy[j*4+1];
				strain_disp_matrix[1][5] = dN_dy[j*4+2];
				strain_disp_matrix[1][7] = dN_dy[j*4+3];

				strain_disp_matrix[2][0] = dN_dy[j*4+0];
				strain_disp_matrix[2][1] = dN_dx[j*4+0];
				strain_disp_matrix[2][2] = dN_dy[j*4+1];
				strain_disp_matrix[2][3] = dN_dx[j*4+1];
				strain_disp_matrix[2][4] = dN_dy[j*4+2];
				strain_disp_matrix[2][5] = dN_dx[j*4+2];
				strain_disp_matrix[2][6] = dN_dy[j*4+3];
				strain_disp_matrix[2][7] = dN_dx[j*4+3];

				// plane stress
				stress_strain_matrix = {
					{1, structure.poisson_ratio, 0},
					{structure.poisson_ratio, 1, 0},
					{0, 0, (1-structure.poisson_ratio)/2}};
				for(k=0;k<3;k++){
					for(l=0;l<3;l++){
						stress_strain_matrix[k][l] *= 1 / (1 - structure.poisson_ratio * structure.poisson_ratio);
					}
				}
				// calculate element stiffness matrix
				calcStiffnessMatrix(stiff_matrix, strain_disp_matrix, stress_strain_matrix);
				calcStiffnessMatrix(stiff_top_matrix, strain_disp_matrix, stress_strain_matrix);
				for(k=0;k<stiff_matrix.size();k++){
					for(l=0;l<stiff_matrix[k].size();l++){
						stiff_matrix[k][l] *= structure.thickness * J[j];
						stiff_top_matrix[k][l] *= structure.thickness * J[j];
						stiff_matrix[k][l] *= structure.youngs_modulus_nodes[node_id[j]];
					}
				}
				// calc strain and stress
				multi_mat_vec(strain_elem, strain_disp_matrix, disp_elem);
				multi_mat_vec(stress_elem, stress_strain_matrix, strain_elem);
				structure.strain_x[elem_id] += strain_elem[0] * J[j];
				structure.strain_y[elem_id] += strain_elem[1] * J[j];
				structure.sheer_strain_xy[elem_id] += strain_elem[2] * J[j];
				structure.stress_x[elem_id] += stress_elem[0] * J[j];
				structure.stress_y[elem_id] += stress_elem[1] * J[j];
				structure.sheer_stress_xy[elem_id] += stress_elem[2] * J[j];
				// calc topology optimization sensitivity
				multi_mat_vec(temp, stiff_top_matrix, disp_elem);
				sens = multi_vec_vec(disp_elem, temp);
				for(k=0;k<4;k++){
					structure.sensitivity[node_id[k]] = sens;
				}
				// set element matrix to adjacency matrix format
				for(k=0;k<4;k++){
					connect_check = 0;
					size_of_column_nonzero = adj_matrix.idx[node_id[k]].size();
					for(l=0;l<4;l++){
						for(m=0;m<size_of_column_nonzero;m++){
							if(adj_matrix.idx[node_id[k]][m] == node_id[l]){
								connect_check++;
								adj_matrix.stiff[node_id[k]*2][m*2]     += stiff_matrix[k*2][l*2];
								adj_matrix.stiff[node_id[k]*2][m*2+1]   += stiff_matrix[k*2][l*2+1];
								adj_matrix.stiff[node_id[k]*2+1][m*2]   += stiff_matrix[k*2+1][l*2];
								adj_matrix.stiff[node_id[k]*2+1][m*2+1] += stiff_matrix[k*2+1][l*2+1];
							}
						}
					}
					if(connect_check != 4){
						std::cerr<<"connectivity error: "<<connect_check<<std::endl;
						exit(1);
					}
				}
			}
		}
#ifdef _OPENMP
		}
#endif
	}
	structure.element_func = f_element_func;
#ifdef MEASURE
	double time_end = elapsedTime();
	std::ofstream ofs(sim_prm.time_output_filename, std::ios::app);
	ofs<<"calcElementMatrix2Dquad() : "<<(time_end - time_start)/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
#ifdef DEBUG
	std::cout<<"Matrix nonzero element num : "<<sim_prm.num_nonzero<<std::endl;
	std::cout<<"Matrix   total element num : "<<(long long)(structure.num_nodes)*structure.num_nodes<<"\n"<<std::endl;
    debugPrintInfo(__func__);
#endif
}

void calcElementMatrix3Dtetra(){
#ifdef MEASURE
	double time_start = elapsedTime();
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
	std::vector<std::vector<double>> mass    = std::vector< std::vector<double> >(4,std::vector<double>(4,0));
	std::vector<std::vector<double>> stiff   = std::vector< std::vector<double> >(4,std::vector<double>(4,0));
	std::vector<std::vector<double>> damping = std::vector< std::vector<double> >(4,std::vector<double>(4,0));
	// fluid local variables for OpenMP parallel
	std::vector<std::vector<double> > f_element_func = structure.element_func;
	// coloring
	int elem_id = 0, color = 0;

	// initialize connectivity format
	allocateAdjMatrix();

	for(color=0;color<(int)structure.colored_elem_id.size();color++){
#ifdef _OPENMP
		#pragma omp parallel
		{
		#pragma omp for firstprivate( \
			j, k, l, volume, connect_check, size_of_column_nonzero, \
			node_id, mat3d, mat4d ,coeff, coeff_x, coeff_y, coeff_z, \
			mass, stiff, damping, \
			elem_id)
#endif 
		for(i=0;i<(int)structure.colored_elem_id[color].size();i++){
			elem_id = structure.colored_elem_id[color][i];
			for(j=0;j<4;j++){
				node_id[j] = structure.element_node_table[elem_id][j];
			}
			// calculate volume of the element
			mat4d = {
				{1.0, 1.0, 1.0, 1.0},
				{structure.x[node_id[0]], structure.x[node_id[1]], structure.x[node_id[2]], structure.x[node_id[3]]},
				{structure.y[node_id[0]], structure.y[node_id[1]], structure.y[node_id[2]], structure.y[node_id[3]]},
				{structure.z[node_id[0]], structure.z[node_id[1]], structure.z[node_id[2]], structure.z[node_id[3]]}
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
					{structure.x[node_id[1]], structure.x[node_id[2]], structure.x[node_id[3]]},
					{structure.y[node_id[1]], structure.y[node_id[2]], structure.y[node_id[3]]},
					{structure.z[node_id[1]], structure.z[node_id[2]], structure.z[node_id[3]]}
				};
				if(j%2==0){
					coeff[j] =   calc_det_3d(mat3d) / 6.0 / volume;
				}else{
					coeff[j] = - calc_det_3d(mat3d) / 6.0 / volume;
				}
				// b_id
				mat3d = {
					{1.0, 1.0, 1.0},
					{structure.y[node_id[1]], structure.y[node_id[2]], structure.y[node_id[3]]},
					{structure.z[node_id[1]], structure.z[node_id[2]], structure.z[node_id[3]]}
				};
				if(j%2==0){
					coeff_x[j] = - calc_det_3d(mat3d) / 6.0 / volume;
				}else{
					coeff_x[j] =   calc_det_3d(mat3d) / 6.0 / volume;
				}
				// c_id
				mat3d = {
					{1.0, 1.0, 1.0},
					{structure.x[node_id[1]], structure.x[node_id[2]], structure.x[node_id[3]]},
					{structure.z[node_id[1]], structure.z[node_id[2]], structure.z[node_id[3]]}
				};
				if(j%2==0){
					coeff_y[j] =   calc_det_3d(mat3d) / 6.0 / volume;
				}else{
					coeff_y[j] = - calc_det_3d(mat3d) / 6.0 / volume;
				}
				// d_id
				mat3d = {
					{1.0, 1.0, 1.0},
					{structure.x[node_id[1]], structure.x[node_id[2]], structure.x[node_id[3]]},
					{structure.y[node_id[1]], structure.y[node_id[2]], structure.y[node_id[3]]}
				};
				if(j%2==0){
					coeff_z[j] = - calc_det_3d(mat3d) / 6.0 / volume;
				}else{
					coeff_z[j] =   calc_det_3d(mat3d) / 6.0 / volume;
				}
			}
			for(j=0;j<4;j++){
				// reset node_id
				node_id[j] = structure.element_node_table[elem_id][j];
			}

			for(j=0;j<4;j++){
				for(k=0;k<4;k++){
					mass[j][k] = (j==k)? volume/10.0 : volume/20.0;
				}
			}
			// set element matrix to adjacency matrix format
			for(j=0;j<4;j++){
				connect_check = 0;
				size_of_column_nonzero = adj_matrix.idx[node_id[j]].size();
				for(k=0;k<4;k++){
					for(l=0;l<size_of_column_nonzero;l++){
						if(adj_matrix.idx[node_id[j]][l] == node_id[k]){
							connect_check++;
							adj_matrix.mass[node_id[j]][l]    += mass[j][k];
							adj_matrix.stiff[node_id[j]][l]   += stiff[j][k];
							adj_matrix.damping[node_id[j]][l] += damping[j][k];
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
	structure.element_func = f_element_func;
#ifdef MEASURE
	double time_end = elapsedTime();
	std::ofstream ofs(sim_prm.time_output_filename, std::ios::app);
	ofs<<"calcElementMatrix3Dtetra(): "<<(time_end - time_start)/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void calcBoundaryShapeFunction(){
	int i = 0, j = 0, node_id;
	std::vector<double> pos_x, pos_y, pos_z;
	double dx1, dy1, dz1, dx2, dy2, dz2, dd1, dd2, d1d2, area = 0;

	// fluid local variables for OpenMP parallel
	for(i=0;i<structure.num_elements;i++){
		// check if element i include boundary surface
		for(j=0;j<4;j++){
			node_id = structure.element_node_table[i][j];
			if(structure.is_boundary[node_id]){
				pos_x.push_back(structure.x[node_id]);
				pos_y.push_back(structure.y[node_id]);
				pos_z.push_back(structure.z[node_id]);
			}
		}
		if(sim_prm.dim == 2 && pos_x.size() == 2){
			dx1 = pos_x[1] - pos_x[0];
			dy1 = pos_y[1] - pos_y[0];
			dz1 = pos_z[1] - pos_z[0];
			area = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1); // area means length in 2D
			for(j=0;j<4;j++){
				node_id = structure.element_node_table[i][j];
				if(structure.is_boundary[node_id]){
					structure.boundary_shape_function[node_id] += area / 2.0;
				}
			}
		}else if(sim_prm.dim == 3 && pos_x.size() == 3){
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
				node_id = structure.element_node_table[i][j];
				if(structure.is_boundary[node_id]){
					structure.boundary_shape_function[node_id] += area / 3.0;
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
void coloringElements(){
	std::vector<std::vector<int>> node2elem = std::vector<std::vector<int>>(structure.num_nodes);
	std::vector<std::vector<int>> elem2node = std::vector<std::vector<int>>(structure.num_elements);
	std::vector<std::vector<int>> elem2elem = std::vector<std::vector<int>>(structure.num_elements);
	std::vector<std::pair<int, int>> deg_size = std::vector<std::pair<int, int>>(structure.num_elements);
	std::vector<int> color_elem = std::vector<int>(structure.num_elements, -1);
	int i, j, k, node_id, elem_id1, elem_id2, max_color = 0;

	// create node-element table
	for(i=0;i<structure.num_elements;i++){
		for(j=0;j<sim_prm.num_polygon_corner;j++){
			node_id = structure.element_node_table[i][j];
			node2elem[node_id].push_back(i);
		}
	}
	// create elemet-node table
	for(i=0;i<structure.num_nodes;i++){
		for(j=0;j<(int)node2elem[i].size();j++){
			elem_id1 = node2elem[i][j];
			elem2node[elem_id1].push_back(i);
		}
	}
	// create element-element table (graph)
	for(i=0;i<structure.num_elements;i++){
		std::vector<bool> is_registered = std::vector<bool>(structure.num_elements, false);
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
	for(i=0;i<structure.num_elements;i++){
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
	structure.colored_elem_id = std::vector<std::vector<int>>(max_color+1);
	for(i=0;i<structure.num_elements;i++){
		int color = color_elem[i];
		structure.colored_elem_id[color].push_back(i);
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
    std::cout<<"Number of colors: "<<max_color<<std::endl;
	for(i=0;i<max_color;i++){
		std::cout<<"Color id: "<<i<<" "<<"Number of nodes: "<<structure.colored_elem_id[i].size()<<std::endl;
	}
#endif
}