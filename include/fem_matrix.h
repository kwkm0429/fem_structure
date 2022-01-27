#pragma once

#include <vector>

void calcJacobian(int, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
void calcElementMatrix2Dquad(void);
void calcElementMatrix3Dtetra(void);
void calcBoundaryShapeFunction(void);
void coloringElements(void);