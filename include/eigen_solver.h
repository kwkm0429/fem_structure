#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Eigenvalues>
#include <vector>

#include "parameter.h"
#include "topopt.h"

typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int64_t> SpMat;
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;

struct StructureMatrix{
	SpMat stiff, stiff_geo;
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
void initSparseMatrix(Sim& sim, Str& str);
void freeSparseMatrix();
void setSparseMatrix(Sim& sim, Str& str, AdjMatrix& adj_mat);
// boundary condition
void setBoundaryCondition2D(SpMat&, Vector&, ParameterID, Sim& sim, Str& str);
// solve linear equation
void solveLinearEquation2D(Sim& sim, Str& str);
// linear solver
bool eigenSolver(SpMat& A, Vector& x, Vector& b, Sim& sim);
bool eigenLU(SpMat& A, Vector& x, Vector& b);
bool eigenBiCGSTAB(SpMat& A, Vector& x, Vector& b);
// buckling solver
void solveBuckling2D(Sim& sim, Str& str);
bool EigenValueSolver(SpMat& A, SpMat& B, Vector& v, double& lambda);
// topology optimization
void calcCompliance(double& compliance);
void calcSensitivity(TopOpt& top, Str& str);