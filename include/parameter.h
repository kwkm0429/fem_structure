#pragma once

#include <vector>
#include <string>

// parameters
struct Sim{
	const std::string time_output_filename = "data-output/time.log";
	const std::string input_data_dirname   = "data-input/";
	const std::string output_data_dirname  = "data-output/";
	const std::string node_filename        = "node.dat";
	const std::string elem_filename        = "elem.dat";
	const std::string bc_filename          = "bc.dat";
	const std::string params_filename      = "sim.prm";

	double dt;
	int max_time_step;
	int cur_time_step = 0;
	int output_interval;
	double gravity_x;
	double gravity_y;
	double gravity_z;
	int eq_solver_opt;
	int num_polygon_corner;
	int num_nonzero;
	int dim;
	int id;
};

struct Str{
	int num_nodes;
	int num_elements;
	double density;
	double youngs_modulus;
	double poisson_ratio;
	double thickness;

	std::vector< std::vector<int> > element_node_table;
	std::vector<double> x, y, z;
	std::vector<double> disp_all, disp_x, disp_y, disp_z;
	// stress and strain
	std::vector<double> stress_x, stress_y, stress_z;
	std::vector<double> sheer_stress_xy, sheer_stress_yz, sheer_stress_zx;
	std::vector<double> strain_x, strain_y, strain_z;
	std::vector<double> sheer_strain_xy, sheer_strain_yz, sheer_strain_zx;
	// boundary conditions
	std::vector<bool> is_boundary;
	std::vector<double> boundary_shape_function;
	std::vector<bool> is_dirichlet_dx, is_dirichlet_dy, is_dirichlet_dz;
	std::vector<double> dirichlet_dx, dirichlet_dy, dirichlet_dz;
	std::vector<double> force_x, force_y, force_z;
	// basis function
	std::vector< std::vector<double> > element_func;
	// coloring
	std::vector<std::vector<int>> colored_elem_id;
	// young's modulus for nodes
	std::vector<double> youngs_modulus_nodes;
	std::vector<double> sensitivity;
};

struct AdjMatrix{
	std::vector< std::vector<int> > idx;
	std::vector< std::vector<double> > stiff, stiff_geo;
};