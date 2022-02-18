#include <vector>
#include <iostream>
#include <fstream>

#include "eigen_solver.h"
#include "time_measure.h"
#include "debug.h"
#include "topopt.h"

static StructureMatrix s_matrix;
static StructureVector s_vector;

void initSparseMatrix(){
#ifdef MEASURE
	double time_start = elapsedTime();
#endif
    // initialize matrices
	s_matrix.mass.resize(structure.num_nodes*sim_prm.dim, structure.num_nodes*sim_prm.dim);
	s_matrix.stiff.resize(structure.num_nodes*sim_prm.dim, structure.num_nodes*sim_prm.dim);
	s_matrix.damping.resize(structure.num_nodes*sim_prm.dim, structure.num_nodes*sim_prm.dim);
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
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void freeSparseMmatrix(){
	// resize
	s_matrix.mass.resize(0, 0);
	s_matrix.stiff.resize(0, 0);
	s_matrix.damping.resize(0, 0);
	// free memory
	s_matrix.mass.data().squeeze();
	s_matrix.stiff.data().squeeze();
	s_matrix.damping.data().squeeze();
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void setSparseMatrix(){
#ifdef MEASURE
	double time_start = elapsedTime();
#endif
	int i,j, nodeIdx;

	for(i=0;i<structure.num_nodes;i++){
		int size_of_column_nonzero = adj_matrix.idx[i].size();
		for(j=0;j<size_of_column_nonzero;j++){
			nodeIdx = adj_matrix.idx[i][j];
			s_matrix.stiff.insert(i*2, nodeIdx*2) = adj_matrix.stiff[i*2][j*2];
			s_matrix.stiff.insert(i*2, nodeIdx*2+1) = adj_matrix.stiff[i*2][j*2+1];
			s_matrix.stiff.insert(i*2+1, nodeIdx*2) = adj_matrix.stiff[i*2+1][j*2];
			s_matrix.stiff.insert(i*2+1, nodeIdx*2+1) = adj_matrix.stiff[i*2+1][j*2+1];
		}
	}
#ifdef MEASURE
	double time_end = elapsedTime();
	std::ofstream ofs(sim_prm.time_output_filename, std::ios::app);
	ofs<<"setSparseMatrix(): "<<(time_end - time_start)/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void setBoundaryCondition2D(SpMat& A, Vector& b){
#ifdef MEASURE
	double time_start = elapsedTime();
#endif
	int i, j;

	for(i=0;i<structure.num_nodes;i++){
		// Neumann Boundary Condition
		b(i*2,0)   += structure.force_x[i];
		b(i*2+1,0) += structure.force_y[i];
		// Dirichlet Boundary Condition
		if(structure.is_dirichlet_dx[i])b(i*2,0)=structure.dirichlet_dx[i];
		if(structure.is_dirichlet_dy[i])b(i*2+1,0)=structure.dirichlet_dy[i];
	}
	for(i=0;i<A.outerSize();++i){
		for(SpMat::InnerIterator it(A, i); it; ++it){
			// Dirichlet Condition
			if(it.col()%2==0){
				if(structure.is_dirichlet_dx[it.col()/2]){
					b(i,0)-=structure.dirichlet_dx[it.col()/2]*it.value();
					if(it.col()==it.row())it.valueRef() = 1;
					else it.valueRef() = 0;
				}
			}else{
				if(structure.is_dirichlet_dy[(it.col()-1)/2]){
					b(i,0)-=structure.dirichlet_dy[(it.col()-1)/2]*it.value();
					if(it.col()==it.row())it.valueRef() = 1;
					else it.valueRef() = 0;
				}
			}
			if(it.row()%2==0){
				if(structure.is_dirichlet_dx[it.row()/2]){
					if(i == it.col()){
						it.valueRef() = 1;
					}else{
						it.valueRef() = 0;
					}
				}
			}else{
				if(structure.is_dirichlet_dy[(it.row()-1)/2]){
					if(i == it.col()){
						it.valueRef() = 1;
					}else{
						it.valueRef() = 0;
					}
				}
			}
		}
	}
#ifdef MEASURE
	double time_end = elapsedTime();
	std::ofstream ofs(sim_prm.time_output_filename, std::ios::app);
	ofs<<"setBoundaryCondition(): "<<(time_end - time_start)/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void solveLinearEquation2D(){
	int i;
	bool is_solved = true;
	Vector rhs = Vector::Zero(structure.num_nodes*2);
	s_vector.disp = Vector::Zero(structure.num_nodes*2);

	// set boundary condition
	setBoundaryCondition2D(s_matrix.stiff, rhs);

	// check matrix and rhs
	//std::cout<<"check stiffness matrix"<<std::endl;
	//std::cout<<s_matrix.stiff<<std::endl;
	//std::cout<<"check rhs vector"<<std::endl;
	//std::cout<<rhs<<std::endl;

	// solve
	is_solved = eigenSolver(s_matrix.stiff, s_vector.disp, rhs);

	for(i=0;i<structure.num_nodes*2;i++){
		structure.disp_all[i] = s_vector.disp(i,0);
		if(i%2==0)structure.disp_x[i/2] = s_vector.disp(i,0);
		else structure.disp_y[(i-1)/2] = s_vector.disp(i,0);
	}

	if(is_solved){
		std::cout<<"Succeeded solveStaticAnalysis\n"<<std::endl;
	}else{
		std::cout<<"Failed to solveStaticAnalysis\n"<<std::endl;
		exit(1);
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

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
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
	return ret;
}

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
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
	return true;
}

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
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
	return true;
}

// Topology Optimization
void calcCompliance(double& compliance){
	compliance = s_vector.disp.transpose() * s_matrix.stiff * s_vector.disp;
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void calcSensitivity(TopOptParameter& top){
	int i,j,node_id;
	Vector stiff_node = Vector::Zero(structure.num_nodes*2);
	Vector sens = Vector::Zero(structure.num_nodes*2);

	for(i=0;i<structure.num_nodes;i++){
		int size_of_column_nonzero = adj_matrix.idx[i].size();
		for(j=0;j<size_of_column_nonzero;j++){
			node_id = adj_matrix.idx[i][j];
			stiff_node(i*2,0)         += adj_matrix.stiff[i*2][j*2];
			//stiff_node(node_id*2,0)   += adj_matrix.stiff[i*2][j*2];
			stiff_node(i*2,0)         += adj_matrix.stiff[i*2][j*2+1];
			//stiff_node(node_id*2+1,0) += adj_matrix.stiff[i*2][j*2+1];
			stiff_node(i*2+1,0)       += adj_matrix.stiff[i*2+1][j*2];
			//stiff_node(node_id*2,0)   += adj_matrix.stiff[i*2+1][j*2];
			stiff_node(i*2+1,0)       += adj_matrix.stiff[i*2+1][j*2+1];
			//stiff_node(node_id*2+1,0) += adj_matrix.stiff[i*2+1][j*2+1];
		}
	}
	sens = s_vector.disp.transpose() * stiff_node * s_vector.disp;
	for(i=0;i<(int)(top.sens.size());i++){
		top.sens[i] = - sens(i,0) * top.pow*(top.E0-top.Emin)*std::pow(top.rho[i],top.pow-1);
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}