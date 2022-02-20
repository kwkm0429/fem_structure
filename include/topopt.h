#pragma once

#include <string>
#include <vector>

#include "parameter.h"

struct TopOpt{
	const std::string params_filename = "top.prm";
	
	std::vector<double> rho;
	std::vector<double> rho_new;
	std::vector<double> sens;
	// input parameter
	double vol_max;
	double rho_init;
	double E0, Emin;
	int pow;
	int itr_max;
	// other parameter
	double vol_sum;
	double vol_frac;
	double comp;
	double comp_prev;
	double comp_conv;
};

void readTopOptDataFile(TopOpt& top, Sim& sim);
void initTopOpt(TopOpt& top, Str& str);
void updateYoungsModulus(TopOpt& top, Str& str);
void calcVolume(TopOpt& top, Str& str);
void optOCmethod(TopOpt& top);
void exeTopOpt(TopOpt& top, Sim& sim, Str& str, AdjMatrix& adj_mat);
