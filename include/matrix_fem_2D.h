#pragma once

#include <vector>

#include "parameter.h"

void calcJacobian(int, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
void calcStiffnessMatrix2D(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
void calcDMatrix2DPlaneStress(std::vector<std::vector<double>>& D, Str& str);
void calcBMatrix2D(std::vector<std::vector<double>>& B, std::vector<double>& dN_dx, std::vector<double>& dN_dy, int j);
void calcDiffMatrix2D(std::vector<std::vector<double>>& Diff, std::vector<double>& dN_dx, std::vector<double>& dN_dy, int j);
void calcStressMatrix2D(std::vector<std::vector<double>>& S, double sigma_x, double sigma_y, double sigma_xy);
void calcStiffGeoMatrix2D(std::vector<std::vector<double>>& Kg, std::vector<std::vector<double>>& Diff, std::vector<std::vector<double>>& S);
void calcMassMatrix2D(std::vector<std::vector<double>>& M, std::vector<double>& N, int j, Str& str);
void calcStressStrain(Sim& sim, Str& str, AdjMatrix& adj_mat);
void calcElementMatrix2Dquad(Sim&, Str& str, AdjMatrix& adj_mat);
void calcBoundaryShapeFunction(Sim&, Str& str);
void coloringElements(Sim&, Str& str);