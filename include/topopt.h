#pragma once

#include <string>
#include <vector>

struct TopOptParameter{
	const std::string params_filename = "top.prm";
	
	std::vector<double> rho;
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

void readTopOptDataFile(void);
void initTopOpt(void);
void updateYoungsModulus(void);
void calcVolume(void);
void optOCmethod(void);
void exeTopOpt(void);
