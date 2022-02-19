#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <vector>

#include "parameter.h"
#include "topopt.h"

typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int64_t> SpMat;
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;

struct StructureMatrix{
	SpMat mass, stiff, damping;
};

struct StructureVector{
	Vector disp;
};

enum ParameterID{
	DX = 0,
	DY = 1,
	DZ = 2,
};

// manage sparse matrix
void initSparseMatrix(void);
void freeSparseMatrix(void);
void setSparseMatrix(void);
// boundary condition
void setBoundaryCondition2D(SpMat&, Vector&, ParameterID);
// solve linear equation
void solveLinearEquation2D(void);
// linear solver
bool eigenSolver(SpMat& A, Vector& x, Vector& b);
bool eigenLU(SpMat& A, Vector& x, Vector& b);
bool eigenBiCGSTAB(SpMat& A, Vector& x, Vector& b);
// topology optimization
void calcCompliance(double& compliance);
void calcSensitivity(TopOptParameter& top);