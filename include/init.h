#pragma once

#include <string>
#include <vector>

void readNodeDataFile(void);
void readElemDataFile(void);
void readBoundaryDataFile(void);
void readParameterDataFile(void);
void initField(void);
void initAdjMatrix(void);
void allocateAdjMatrix(void);