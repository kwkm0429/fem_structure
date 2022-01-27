#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <vector>

#include "parameter.h"

typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int64_t> SpMat;
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;

struct StructureMatrix{
	SpMat mass, stiff, damping;
};

enum ParameterID{
	VX = 0,
	VY = 1,
	VZ = 2,
	PRES = 3
};

// manage sparse matrix
void initSparseMatrix(void);
void freeSparseMmatrix(void);
void setSparseMatrix(void);
// boundary condition
void setBoundaryCondition(SpMat&, Vector&, ParameterID);
// linear solver
bool eigenSolver(SpMat& A, Vector& x, Vector& b);
bool eigenLU(SpMat& A, Vector& x, Vector& b);
bool eigenBiCGSTAB(SpMat& A, Vector& x, Vector& b);