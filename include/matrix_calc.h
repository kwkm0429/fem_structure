#pragma once

#include <vector>

double multi_vec_vec(const std::vector<double>&, const std::vector<double>&,int );
void multi_mat_vec(std::vector<double>&, const std::vector< std::vector<double> >&, const std::vector<double>&,int);
void multi_mat_mat(std::vector< std::vector<double> >&, const std::vector< std::vector<double> >&, const std::vector< std::vector<double> >&);
void transpose_mat(std::vector< std::vector<double> >&, const std::vector< std::vector<double> >&);
bool inv_mat_2d(std::vector< std::vector<double> >&, const std::vector< std::vector<double> >&);
bool inv_mat_3d(std::vector< std::vector<double> >&, const std::vector< std::vector<double> >&);
bool inv_mat_4d(std::vector< std::vector<double> >&, const std::vector< std::vector<double> >&);
double calc_trace(const std::vector< std::vector<double> >&, int);
double calc_det_2d(const std::vector< std::vector<double> >&);
double calc_det_3d(const std::vector< std::vector<double> >&);
double calc_det_4d(const std::vector< std::vector<double> > &);