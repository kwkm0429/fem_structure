#pragma once

#include "parameter.h"

void updatePosition(Str& str);
void readInputFiles(Sim& sim, Str& str);
void initStructureStatus(Sim& sim, Str& str, AdjMatrix& adj_mat);
void exePostProcess(Sim& sim, Str& str, AdjMatrix& adj_mat);
void exeStaticAnalysis(Sim& sim, Str& str, AdjMatrix& adj_mat);
void exeBucklingAnalysis(Sim& sim, Str& str, AdjMatrix& adj_mat);
void exeModalAnalysis(Sim& sim, Str& str, AdjMatrix& adj_mat);