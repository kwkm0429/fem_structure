#include <vector>
#include <iostream>
#include <fstream>

#include "eigen_solver.h"
#include "time_measure.h"
#include "debug.h"
#include "topopt.h"

static StructureMatrix s_mat;
static StructureVector s_vec;

void initSparseMatrix(Sim& sim, Str& str){
#ifdef MEASURE
	double time_start = elapsedTime();
#endif
    // initialize matrices
	s_mat.stiff.resize(str.num_nodes*sim.dim, str.num_nodes*sim.dim);
	s_mat.stiff_geo.resize(str.num_nodes*sim.dim, str.num_nodes*sim.dim);
	// reserve memory
	s_mat.stiff.reserve(sim.num_nonzero);
	s_mat.stiff_geo.reserve(sim.num_nonzero);
#ifdef MEASURE
	double time_end = elapsedTime();
	std::ofstream ofs(sim.time_output_filename, std::ios::app);
	ofs<<"initSparseMatrix(): "<<(time_end - time_start)/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void freeSparseMatrix(){
	// resize
	s_mat.stiff.resize(0, 0);
	s_mat.stiff_geo.resize(0, 0);
	// free memory
	s_mat.stiff.data().squeeze();
	s_mat.stiff_geo.data().squeeze();
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void setSparseMatrix(Sim& sim, Str& str, AdjMatrix& adj_mat){
#ifdef MEASURE
	double time_start = elapsedTime();
#endif
	int i,j, nodeIdx;
	for(i=0;i<str.num_nodes*2;i++){
		int size_of_column_nonzero = adj_mat.idx[(int)(i/2)].size();
		for(j=0;j<size_of_column_nonzero;j++){
			nodeIdx = adj_mat.idx[(int)(i/2)][j];
			s_mat.stiff.insert(i, nodeIdx*2) = adj_mat.stiff[i][j*2];
			s_mat.stiff.insert(i, nodeIdx*2+1) = adj_mat.stiff[i][j*2+1];
			s_mat.stiff_geo.insert(i, nodeIdx*2) = adj_mat.stiff_geo[i][j*2];
			s_mat.stiff_geo.insert(i, nodeIdx*2+1) = adj_mat.stiff_geo[i][j*2+1];
		}
	}
#ifdef MEASURE
	double time_end = elapsedTime();
	std::ofstream ofs(sim.time_output_filename, std::ios::app);
	ofs<<"setSparseMatrix(): "<<(time_end - time_start)/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void setBoundaryCondition2D(SpMat& A, Vector& b, Sim& sim, Str& str){
#ifdef MEASURE
	double time_start = elapsedTime();
#endif
	int i, j;

	for(i=0;i<str.num_nodes;i++){
		// Neumann Boundary Condition
		b(i*2,0)   += str.force_x[i];
		b(i*2+1,0) += str.force_y[i];
		// Dirichlet Boundary Condition
		if(str.is_dirichlet_dx[i])b(i*2,0)=str.dirichlet_dx[i];
		if(str.is_dirichlet_dy[i])b(i*2+1,0)=str.dirichlet_dy[i];
	}
	for(i=0;i<A.outerSize();++i){
		for(SpMat::InnerIterator it(A, i); it; ++it){
			// Dirichlet Condition
			if(it.col()%2==0){
				if(str.is_dirichlet_dx[it.col()/2]){
					b(i,0)-=str.dirichlet_dx[it.col()/2]*it.value();
					if(it.col()==it.row())it.valueRef() = 1;
					else it.valueRef() = 0;
				}
			}else{
				if(str.is_dirichlet_dy[(it.col()-1)/2]){
					b(i,0)-=str.dirichlet_dy[(it.col()-1)/2]*it.value();
					if(it.col()==it.row())it.valueRef() = 1;
					else it.valueRef() = 0;
				}
			}
			if(it.row()%2==0){
				if(str.is_dirichlet_dx[it.row()/2]){
					if(i == it.col()){
						it.valueRef() = 1;
					}else{
						it.valueRef() = 0;
					}
				}
			}else{
				if(str.is_dirichlet_dy[(it.row()-1)/2]){
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
	std::ofstream ofs(sim.time_output_filename, std::ios::app);
	ofs<<"setBoundaryCondition(): "<<(time_end - time_start)/1000<<" [s]"<<std::endl;
	ofs.close();
#endif
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void solveLinearEquation2D(Sim& sim, Str& str){
	int i;
	bool is_solved = true;
	Vector rhs = Vector::Zero(str.num_nodes*2);
	s_vec.disp = Vector::Zero(str.num_nodes*2);

	// set boundary condition
	setBoundaryCondition2D(s_mat.stiff, rhs, sim, str);

	// check matrix and rhs
	//std::cout<<"check stiffness matrix"<<std::endl;
	//std::cout<<s_mat.stiff<<std::endl;
	//std::cout<<"check rhs vector"<<std::endl;
	//std::cout<<rhs<<std::endl;

	// solve
	is_solved = eigenSolver(s_mat.stiff, s_vec.disp, rhs, sim);

	for(i=0;i<str.num_nodes*2;i++){
		str.disp_all[i] = s_vec.disp(i,0);
		if(i%2==0)str.disp_x[i/2] = s_vec.disp(i,0);
		else str.disp_y[(i-1)/2] = s_vec.disp(i,0);
	}

	if(!is_solved){
		std::cout<<"Failed "<<__func__<<std::endl;
		exit(1);
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

bool eigenSolver(SpMat& A, Vector& x, Vector& b, Sim& sim){
#ifdef MEASURE
	double time_start = elapsedTime();
#endif 
#ifdef _OPENMP
	Eigen::initParallel();
#endif
	bool ret;
	switch(sim.eq_solver_opt){
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
	std::ofstream ofs(sim.time_output_filename, std::ios::app);
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

void solveBuckling2D(Sim& sim, Str& str){
	int i;
	bool is_solved = true;
	Vector v = Vector::Zero(str.num_nodes*2);

	// solve (KL+\lambda KG)\phi = 0
	//std::cout<<"check geometric stiffness matrix"<<std::endl;
	//std::cout<<s_mat.stiff_geo<<std::endl;
	is_solved = EigenValueSolver(s_mat.stiff, s_mat.stiff_geo, v, str.buckling_coeff);
	for(i=0;i<str.num_nodes;i++){
		str.buckling_x[i] = v[i*2];
		str.buckling_y[i] = v[i*2+1];
	}
	if(!is_solved){
		std::cout<<"Failed "<<__func__<<std::endl;
		exit(1);
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

bool EigenValueSolver(SpMat& A, SpMat& B, Vector& v, double& lambda){
	int i, id, size;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A,B);
	//std::cout << "The eigenvalues of the pencil (A,B) are:" << std::endl << es.eigenvalues() << std::endl;
	//std::cout << "The matrix of eigenvectors, V, is:" << std::endl << es.eigenvectors() << std::endl << std::endl;
	
	size = es.eigenvalues().size();
	for(i=0;i<size;i++){
		if(es.eigenvalues()[i]>1e-6){
			lambda=es.eigenvalues()[i];
			id = i;
			break;
		}
	}
	std::cout << "Consider the first eigenvalue, lambda = " << lambda << std::endl;
	v = es.eigenvectors().col(id);
	//std::cout << "If v is the corresponding eigenvector, then A * v = " << std::endl << A * v << std::endl;
	//std::cout << "... and lambda * B * v = " << std::endl << lambda * B * v << std::endl << std::endl;
	if(es.info()!=Eigen::Success) {
		// solving failed
		std::cout<<"Failed "<<__func__<<std::endl;
		return false;
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
	return true;
}

// Topology Optimization
void calcCompliance(double& compliance){
	compliance = s_vec.disp.transpose() * s_mat.stiff * s_vec.disp;
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void calcSensitivity(TopOpt& top, Str& str){
	int i,j,node_id;
	for(i=0;i<(int)(top.sens.size());i++){
		top.sens[i] = - str.sensitivity[i]*top.pow*(top.E0-top.Emin)*std::pow(top.rho[i],top.pow-1);
		//std::cout<<str.sensitivity[i]<<" "<<top.sens[i]<<" "<<str.youngs_modulus_nodes[i]<<" "<<top.rho[i]<<std::endl;
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}