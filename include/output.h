#pragma once

#include "parameter.h"
#include "topopt.h"

void outputVtkFile(int, Sim& sim, Str& str);
void outputDispVtkFile(int, Sim& sim, Str& str);
void outputStrainVtkFile(int, Sim& sim, Str& str);
void outputStressVtkFile(int, Sim& sim, Str& str);
void outputBucklingVtkFile(int, Sim& sim, Str& str);
void outputEigenModeVtkFile(int number, Sim& sim, Str& str);
void outputDensityVtkFile(int, TopOpt&, Sim& sim, Str& str);
void outputParameterDataFile(Sim& sim);
void outputTopOptDataFile(TopOpt&, Sim& sim);