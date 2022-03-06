#pragma once

#include <vector>

#include "parameter.h"

void calcStiffnessMatrix3D(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
void calcDMatrix3DPlaneStress(std::vector<std::vector<double>>& D, Str& str);
void calcBMatrix3D(std::vector<std::vector<double>>& B, std::vector<double>& dN_dx, std::vector<double>& dN_dy, int j);
void calcDiffMatrix3D(std::vector<std::vector<double>>& Diff, std::vector<double>& dN_dx, std::vector<double>& dN_dy, int j);
void calcStressMatrix3D(std::vector<std::vector<double>>& S, double sigma_x, double sigma_y, double sigma_xy);
void calcStiffGeoMatrix3D(std::vector<std::vector<double>>& Kg, std::vector<std::vector<double>>& Diff, std::vector<std::vector<double>>& S);
void calcMassMatrix3D(std::vector<std::vector<double>>& M, std::vector<double>& N, int j, Str& str);
void calcElementMatrix3Dtetra(Sim&, Str& str, AdjMatrix& adj_mat);