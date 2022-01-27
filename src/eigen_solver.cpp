/*
* @Author: Kosuke Kawakami
* @Date:   2019-08-06 16:08:01
* @Last Modified by:   Kosuke
* @Last Modified time: 2020-08-23 11:20:18
*/
#include <vector>
#include <iostream>
#include <fstream>

#include "eigen_solver.h"
#include "time_measure.h"

static StructureMatrix s_matrix;

/**
 * @brief      Initialize sparse matrix (s_matrix) and set memory.
 * 
 */
void initSparseMatrix(){
#ifdef MEASURE
	double time_start = elapsedTime();
#endif
    // initialize matrices
	s_matrix.mass.resize(structure.num_nodes, structure.num_nodes);
	s_matrix.stiff.resize(structure.num_nodes, structure.num_nodes);
	s_matrix.damping.resize(structure.num_nodes, structure.num_nodes);
	// reserve memory
	s_matrix.mass.reserve(sim_prm.num_nonzero);
	s_matrix.stiff.reserve(sim_prm.num_nonzero);
	s_matrix.damping.reserve(sim_prm.num_nonzero);
#ifdef MEASURE
	double time_end = elapsedTime();
	std::ofstream ofs(sim_prm.time_output_filename, std::ios::app);
	ofs<<"initSparseMatrix(): "<<(time_end - time_start)/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
}

/**
 * @brief      Free memory of sparse matrix (s_matrix)
 */
void freeSparseMmatrix(){
	// resize
	s_matrix.mass.resize(0, 0);
	s_matrix.stiff.resize(0, 0);
	s_matrix.damping.resize(0, 0);
	// free memory
	s_matrix.mass.data().squeeze();
	s_matrix.stiff.data().squeeze();
	s_matrix.damping.data().squeeze();
}

/**
 * @brief      Set values of elements in sparse matrix
 *
 * @param      idx        index of nodes
 * @param      mass       mass matrix
 * @param      diffuse    diffuse matrix
 * @param      pressure   pressure matrix
 * @param      conti_x    continuity matrix in x
 * @param      conti_y    continuity matrix in y 
 * @param      conti_z    continuity matrix in z 
 * @param      mass_stab  mass stabilization matrix
 * @param      adv        advection matrix
 * @param      adv_stab   advection stabilization matrix
 */
void setSparseMatrix(){
#ifdef MEASURE
	double time_start = elapsedTime();
#endif
	int i,j, nodeIdx;

	for(i=0;i<structure.num_nodes;i++){
		int size_of_column_nonzero = adj_matrix.idx[i].size();
		for(j=0;j<size_of_column_nonzero;j++){
			nodeIdx = adj_matrix.idx[i][j];
			s_matrix.mass.insert(i, nodeIdx)    = adj_matrix.mass[i][j];
			s_matrix.stiff.insert(i, nodeIdx)   = adj_matrix.stiff[i][j];
			s_matrix.damping.insert(i, nodeIdx) = adj_matrix.damping[i][j];
		}
	}
#ifdef MEASURE
	double time_end = elapsedTime();
	std::ofstream ofs(sim_prm.time_output_filename, std::ios::app);
	ofs<<"setSparseMatrix(): "<<(time_end - time_start)/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
}

/**
 * @brief      Set the boundary condition for velocity in x.
 *
 * @param      A    { Coefficient Matrix }
 * @param      b    { Right hand side vector }
 * @param      param {VX, VY, VZ, PRES}
 */
void setBoundaryCondition(SpMat& A, Vector& b, ParameterID param_id){
#ifdef MEASURE
	double time_start = elapsedTime();
#endif
	int i;

	for(i=0;i<structure.num_nodes;i++){
		if(param_id == VX){
			if(structure.is_dirichlet_vx[i]){
				// Dirichlet Boundary Condition
				b(i,0) = structure.dirichlet_vx[i] * structure.density;
			}else{
				// Neumann Boundary Condition
				b(i,0) += (structure.neumann_vx[i][0] + structure.neumann_vx[i][1] + structure.neumann_vx[i][2]) * structure.visc / structure.density * structure.boundary_shape_function[i];
			}
		}else if(param_id == VY){
			if(structure.is_dirichlet_vy[i]){
				// Dirichlet Boundary Condition
				b(i,0) = structure.dirichlet_vy[i] * structure.density;
			}else{
				// Neumann Boundary Condition
				b(i,0) += (structure.neumann_vy[i][0] + structure.neumann_vy[i][1] + structure.neumann_vy[i][2]) * structure.visc / structure.density * structure.boundary_shape_function[i];
			}
		}else if(param_id == VZ){
			if(structure.is_dirichlet_vz[i]){
				// Dirichlet Boundary Condition
				b(i,0) = structure.dirichlet_vz[i] * structure.density;

			}else{
				// Neumann Boundary Condition
				b(i,0) += (structure.neumann_vz[i][0] + structure.neumann_vz[i][1] + structure.neumann_vz[i][2]) * structure.visc / structure.density * structure.boundary_shape_function[i];
			}
		}else if(param_id == PRES){
			if(structure.is_dirichlet_pressure[i]){
				// Dirichlet Boundary Condition
				b(i,0) = structure.dirichlet_pressure[i];
			}else{
				// Neumann Boundary Condition
				b(i,0) += (structure.neumann_pressure[i][0] + structure.neumann_pressure[i][1] + structure.neumann_pressure[i][2]) * structure.boundary_shape_function[i];
			}
		}else{
			std::cout<<"Error: "<<param_id<<" , ParameterID error in setBoundaryCondition()"<<std::endl;
			exit(1);
		}
	}
	for(i=0;i<A.outerSize();++i){
		for(SpMat::InnerIterator it(A, i); it; ++it){
			if(param_id == VX){
				if(structure.is_dirichlet_vx[it.row()]){
					if(i == it.col()){
						it.valueRef() = 1;
					}else{
						it.valueRef() = 0;
					}
				}
			}else if(param_id == VY){
				if(structure.is_dirichlet_vy[it.row()]){
					if(i == it.col()){
						it.valueRef() = 1;
					}else{
						it.valueRef() = 0;
					}
				}
			}else if(param_id == VZ){
				if(structure.is_dirichlet_vz[it.row()]){
					if(i == it.col()){
						it.valueRef() = 1;
					}else{
						it.valueRef() = 0;
					}
				}
			}else if(param_id == PRES){
				if(structure.is_dirichlet_pressure[it.row()]){
					if(i == it.col()){
						it.valueRef() = 1;
					}else{
						it.valueRef() = 0;
					}
				}
			}else{
				std::cout<<"Error: "<<param_id<<" , ParameteriD error in setBoundaryCondition()"<<std::endl;
				exit(1);
			}
		}
	}
#ifdef MEASURE
	double time_end = elapsedTime();
	std::ofstream ofs(sim_prm.time_output_filename, std::ios::app);
	ofs<<"setBoundaryCondition(): "<<(time_end - time_start)/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
}

/**
 * @brief      Choose linear equation solution algorithm
 *
 * @param      A     { Coefficient Matrix }
 * @param      x     { Unknoen parameter vector }
 * @param      b     { Right hand side vector }
 *
 * @return     { if succeeded in solving or not }
 */
bool eigenSolver(SpMat& A, Vector& x, Vector& b){
#ifdef MEASURE
	double time_start = elapsedTime();
#endif 
#ifdef _OPENMP
	Eigen::initParallel();
#endif
	bool ret;
	switch(sim_prm.eq_solver_opt){
		case 0:
			ret = eigenLU(A, x, b);
			break;
		case 1:
			ret = eigenBiCGSTAB(A, x, b);
			break;
		default:
			std::cerr<<"EQUATION SOLVER OPTION ERROR"<<std::endl;
			ret = false;
	}
#ifdef MEASURE
	double time_end = elapsedTime();
	std::ofstream ofs(sim_prm.time_output_filename, std::ios::app);
	ofs<<"eigenSolver(): "<<(time_end - time_start)/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
	return ret;
}

/**
 * @brief      Solve linear equation with LU
 *
 * @param      A     { Coefficient Matrix }
 * @param      x     { Unknoen parameter vector }
 * @param      b     { Right hand side vector }
 *
 * @return     { if succeeded in solving or not }
 */
bool eigenLU(SpMat& A, Vector& x, Vector& b){
	Eigen::SparseLU<SpMat> solver;
	solver.compute(A);
	if(solver.info()!=Eigen::Success) {
		// decomposition failed
		std::cout<<"Decomposition failed"<<std::endl;
		return false;
	}
	x = solver.solve(b);
	if(solver.info()!=Eigen::Success) {
		// solving failed
		std::cout<<"Failed solver LU"<<std::endl;
		return false;
	}
	return true;
}

/**
 * @brief      { Solve linear equation with BiCGSTAB }
 *
 * @param      A     { Coefficient Matrix }
 * @param      x     { Unknoen parameter vector }
 * @param      b     { Right hand side vector }
 *
 * @return     { if succeeded in solving or not }
 */
bool eigenBiCGSTAB(SpMat& A, Vector& x, Vector& b){
	Eigen::BiCGSTAB<SpMat> solver;
	solver.setTolerance(1e-9);
	//solver.setMaxIterations(10000);
	solver.compute(A);
	if(solver.info()!=Eigen::Success) {
		// decomposition failed
		std::cout<<"Decomposition failed"<<std::endl;
		return false;
	}
	x = solver.solve(b);
	if(solver.info()!=Eigen::Success) {
		// solving failed
		std::cout<<"Failed solver BiCGSTAB"<<std::endl;
		return false;
	}
	//std::cout << "#nbThreads:       " << omp_get_max_threads() <<std::endl;
	std::cout << "#iterations:     " << solver.iterations() << std::endl;
	std::cout << "estimated error: " << solver.error()      << std::endl;
	return true;
}