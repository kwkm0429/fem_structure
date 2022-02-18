#pragma once

#include "topopt.h"

void outputVtkFile(int);
void outputDispVtkFile(int);
void outputStrainVtkFile(int);
void outputStressVtkFile(int);
void outputDensityVtkFile(int, TopOptParameter&);
void outputParameterDataFile(void);
void outputTopOptDataFile(TopOptParameter& top);