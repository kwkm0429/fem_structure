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
	double t_start = elapsedTime();
#endif
    // initialize matrices
	s_mat.stiff.resize(str.num_nodes*sim.dim, str.num_nodes*sim.dim);
	s_mat.stiff_geo.resize(str.num_nodes*sim.dim, str.num_nodes*sim.dim);
	s_mat.mass.resize(str.num_nodes*sim.dim, str.num_nodes*sim.dim);
	// reserve memory
	s_mat.stiff.reserve(sim.num_nonzero);
	s_mat.stiff_geo.reserve(sim.num_nonzero);
	s_mat.mass.reserve(sim.num_nonzero);
#ifdef MEASURE
	double t_end = elapsedTime();
	writeTime(sim.time_output_filename, __func__, t_start, t_end);
#endif
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void freeSparseMatrix(){
	// resize
	s_mat.stiff.resize(0, 0);
	s_mat.stiff_geo.resize(0, 0);
	s_mat.mass.resize(0, 0);
	// free memory
	s_mat.stiff.data().squeeze();
	s_mat.stiff_geo.data().squeeze();
	s_mat.mass.data().squeeze();
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void setSparseMatrix(Sim& sim, Str& str, AdjMatrix& adj_mat){
#ifdef MEASURE
	double t_start = elapsedTime();
#endif
	int i,j, nodeIdx;
	for(i=0;i<str.num_nodes*2;i++){
		int size_of_column_nonzero = adj_mat.idx[(int)(i/2)].size();
		for(j=0;j<size_of_column_nonzero;j++){
			nodeIdx = adj_mat.idx[(int)(i/2)][j];
			s_mat.stiff.insert(i, nodeIdx*2) = adj_mat.stiff[i][j*2];
			s_mat.stiff.insert(i, nodeIdx*2+1) = adj_mat.stiff[i][j*2+1];
			s_mat.stiff_geo.insert(i, nodeIdx*2) = -adj_mat.stiff_geo[i][j*2];
			s_mat.stiff_geo.insert(i, nodeIdx*2+1) = -adj_mat.stiff_geo[i][j*2+1];
			s_mat.mass.insert(i, nodeIdx*2) = adj_mat.mass[i][j*2];
			s_mat.mass.insert(i, nodeIdx*2+1) = adj_mat.mass[i][j*2+1];
		}
	}
#ifdef MEASURE
	double t_end = elapsedTime();
	writeTime(sim.time_output_filename, __func__, t_start, t_end);
#endif
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void setBoundaryCondition2D(SpMat& A, Vector& b, Sim& sim, Str& str){
#ifdef MEASURE
	double t_start = elapsedTime();
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
	double t_end = elapsedTime();
	writeTime(sim.time_output_filename, __func__, t_start, t_end);
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
	double t_start = elapsedTime();
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
	double t_end = elapsedTime();
	writeTime(sim.time_output_filename, __func__, t_start, t_end);
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

void solveEigenMode2D(Sim& sim, Str& str){
	int i, j;
	bool is_solved = true;
	Vector v = Vector::Zero(str.num_nodes*2);
	Matrix phi = Matrix::Zero(str.num_nodes*2, sim.num_mode);

	//std::cout<<"check geometric stiffness matrix"<<std::endl;
	//std::cout<<s_mat.mass<<std::endl;
	//std::cout<<s_mat.stiff<<std::endl;
	setBoundaryCondition2D(s_mat.mass, v, sim ,str);
	setBoundaryCondition2D(s_mat.stiff, v, sim ,str);
	// solve Kv = \lambda Mv
	is_solved = EigenValueSolver(s_mat.stiff, s_mat.mass, phi, str.eigen_value, sim);
	//is_solved = SpectraSolver(s_mat.stiff_geo, s_mat.stiff, v, str.buckling_coeff, sim);

	for(i=0;i<sim.num_mode;i++){
		for(j=0;j<str.num_nodes;j++){
			str.eigen_mode_x[i][j] = phi(j*2, i);
			str.eigen_mode_y[i][j] = phi(j*2+1, i);
		}
	}
	if(!is_solved){
		std::cout<<"Failed "<<__func__<<std::endl;
		exit(1);
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}


void solveBuckling2D(Sim& sim, Str& str){
	int i, j;
	bool is_solved = true;
	Vector v = Vector::Zero(str.num_nodes*2);
	Matrix phi = Matrix::Zero(str.num_nodes*2, sim.num_mode);

	//std::cout<<"check geometric stiffness matrix"<<std::endl;
	//std::cout<<s_mat.stiff<<std::endl;
	//std::cout<<s_mat.stiff_geo<<std::endl;
	setBoundaryCondition2D(s_mat.stiff_geo, v, sim ,str);
	setBoundaryCondition2D(s_mat.stiff, v, sim ,str);
	// solve (KL+\lambda Kg)\phi = 0 -> (Kg+1/\lambda KL)\phi = 0 -> -Kg v = 1/lambda KL v
	//is_solved = EigenValueSolver(s_mat.stiff_geo, s_mat.stiff, phi, str.buckling_coeff, sim);
	is_solved = SpectraSolver(s_mat.stiff_geo, s_mat.stiff, phi, str.buckling_coeff, sim);

	for(i=0;i<sim.num_mode;i++){
		for(j=0;j<str.num_nodes;j++){
			str.buckling_x[i][j] = phi(j*2, i);
			str.buckling_y[i][j] = phi(j*2+1, i);
		}
	}
	if(!is_solved){
		std::cout<<"Failed "<<__func__<<std::endl;
		exit(1);
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

bool EigenValueSolver(SpMat& A, SpMat& B, Matrix& phi, std::vector<double>& lambda, Sim& sim){
#ifdef MEASURE
	double t_start = elapsedTime();
#endif
	int i, j, id, size;
	double norm=0;
	// convert sparse to dense matrix because Eigen does not support sparse matrix for eigenvalue solver.
	Matrix Ad = A;
	Matrix Bd = B;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(Ad, Bd);
	std::cout << "The eigenvalues of the pencil (A,B) are:" << std::endl << es.eigenvalues() << std::endl;
	//std::cout << "The matrix of eigenvectors, V, is:" << std::endl << es.eigenvectors() << std::endl << std::endl;
	
	size = es.eigenvalues().size();
	id = 0;
	for(i=0;id<sim.num_mode;i++){
		if(std::abs(es.eigenvalues()[size-1-i]-1)<1e-7)continue;
		if(id>=sim.num_mode)break;
		lambda[id] = 1.0/es.eigenvalues()[size-1-i];
		phi.col(id) = es.eigenvectors().col(size-1-i);
		std::cout << id<<") Eigen value, lambda = " << lambda[id] << std::endl;
		id++;
	}
	if(es.info()!=Eigen::Success) {
		// solving failed
		std::cout<<"Failed "<<__func__<<std::endl;
		return false;
	}
#ifdef MEASURE
	double t_end = elapsedTime();
	writeTime(sim.time_output_filename, __func__, t_start, t_end);
#endif
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
	return true;
}

/*
bool SpectraSolver(SpMat& A, SpMat& B, Vector& v, std::vector<double>& lambda, Sim& sim){
	int i;
	// convert to Eigen::ColMajor because Spectra and Eigen has error for RowMajor
	//*
	Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> Ae, Be;
	Ae.resize(A.outerSize(), A.outerSize());
	Be.resize(B.outerSize(), B.outerSize());
	Ae.reserve(sim.num_nonzero);
	Be.reserve(sim.num_nonzero);
	for(i=0;i<A.outerSize();++i){
		for(SpMat::InnerIterator it(A, i); it; ++it){
			Ae.insert(it.col(), i) = it.valueRef();
		}
	}
	for(i=0;i<B.outerSize();++i){
		for(SpMat::InnerIterator it(B, i); it; ++it){
			Be.insert(it.col(), i) = it.valueRef();
		}
	}
	std::cout<<"GATE1"<<std::endl;
	using OpType = Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse, Eigen::Upper, Eigen::Upper, Eigen::RowMajor, Eigen::RowMajor, int64_t, int64_t>;
	using BOpType = Spectra::SparseSymMatProd<double, Eigen::Upper, Eigen::RowMajor, int64_t>;
	OpType op(A, B);
	BOpType Bop(A);
	std::cout<<"GATE2"<<std::endl;
	// Construct generalized eigen solver object, seeking three generalized
    // eigenvalues that are closest to and larger than 1.0. This is equivalent to
    // specifying a shift sigma = 1.0 combined with the SortRule::LargestAlge
    // selection rule
	Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::Buckling> geigs(op, Bop, 5, 10, 1.0);
	std::cout<<"GATE3"<<std::endl;
	geigs.init();
	int nconv = geigs.compute(Spectra::SortRule::LargestAlge);
	std::cout<<"GATE4"<<std::endl;
	if(Spectra::CompInfo::Successful != geigs.info()){
		std::cerr << "Failed " << __func__ <<std::endl;
		exit(-1);
	}else{
		Vector evalues = geigs.eigenvalues();
		std::cout<<evalues<<std::endl;
		Matrix evecs = geigs.eigenvectors();
		v = evecs.col(0);
		lambda[0] = 1.0/evalues(0,0);
		std::cout << "Consider the first eigenvalue, lambda = " << lambda[0] << std::endl;
		return true;
	}
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}*/

//*
bool SpectraSolver(SpMat& A, SpMat& B, Matrix& phi, std::vector<double>& lambda, Sim& sim){
#ifdef MEASURE
	double t_start = elapsedTime();
#endif
	int i, j, id, size;
	// convert to Eigen::ColMajor because Spectra and Eigen has error for RowMajor
	//*
	Eigen::SparseMatrix<double> Ae, Be;
	Ae.resize(A.outerSize(), A.outerSize());
	Be.resize(B.outerSize(), B.outerSize());
	Ae.reserve(sim.num_nonzero);
	Be.reserve(sim.num_nonzero);
	for(i=0;i<A.outerSize();++i){
		for(SpMat::InnerIterator it(A, i); it; ++it){
			Ae.insert(it.col(), i) = it.valueRef();
		}
	}
	for(i=0;i<B.outerSize();++i){
		for(SpMat::InnerIterator it(B, i); it; ++it){
			Be.insert(it.col(), i) = it.valueRef();
		}
	}
	Spectra::SparseSymMatProd<double> Aop(Ae);
	Spectra::SparseCholesky<double> Bop(Be);
	Spectra::SymGEigsSolver<Spectra::SparseSymMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEigsMode::Cholesky> geigs(Aop, Bop, sim.num_mode+sim.num_spc, sim.num_mode*2+sim.num_spc);
	geigs.init();
	int nconv = geigs.compute(Spectra::SortRule::LargestAlge);

	size = geigs.eigenvalues().size();
	id = 0;
	std::cout<<geigs.eigenvalues()<<std::endl;
	for(i=sim.num_spc;i<sim.num_spc+sim.num_mode;i++){
		lambda[id] = 1.0/geigs.eigenvalues()[i];
		phi.col(id) = geigs.eigenvectors().col(i);
		std::cout << id<<") Eigen value, lambda = " << lambda[id] << std::endl;
		id++;
	}

	if(Spectra::CompInfo::Successful != geigs.info()){
		std::cerr << "Failed " << __func__ <<std::endl;
		return false;
	}
#ifdef MEASURE
	double t_end = elapsedTime();
	writeTime(sim.time_output_filename, __func__, t_start, t_end);
#endif
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
    return true;
}//*/

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