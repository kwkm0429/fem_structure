#pragma once

#include <vector>

#include "parameter.h"

void calcJacobian(int, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
void calcStiffnessMatrix(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
void calcElementMatrix2Dquad(Sim&, Str& str, AdjMatrix& adj_mat);
void calcElementMatrix3Dtetra(Sim&, Str& str, AdjMatrix& adj_mat);
void calcBoundaryShapeFunction(Sim&, Str& str);
void coloringElements(Sim&, Str& str);