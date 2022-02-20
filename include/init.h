#pragma once

#include <string>
#include <vector>

#include "parameter.h"

void readNodeDataFile(Sim& sim, Str& str);
void readElemDataFile(Sim& sim, Str& str);
void readBoundaryDataFile(Sim& sim, Str& str);
void readParameterDataFile(Sim& sim, Str& str);
void initField(Sim& sim, Str& str);
void initAdjMatrix(Sim& sim, Str& str, AdjMatrix& adj_mat);
void allocateAdjMatrix(Sim& sim, Str& str, AdjMatrix& adj_mat);