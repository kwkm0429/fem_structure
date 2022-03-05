#pragma once

#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>

#define MEASURE

void initializeTime();
double elapsedTime();
void writeTime(std::string file_name, std::string func, double t_start, double t_end);