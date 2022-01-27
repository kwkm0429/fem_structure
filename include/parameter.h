#pragma once

#include <vector>
#include <string>

// parameters
struct SimulationParameter{
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
};

struct StructureParameter{
	int num_nodes;
	int num_elements;
	double visc;
	double density;
	double reynolds;
	double youngs_modulus;
	double poisson_ratio;

	std::vector< std::vector<int> > element_node_table;
	std::vector<double> x, y, z;
	std::vector<double> disp_x, disp_y, disp_z;
	// stress and strain
	std::vector<double> stress_x, stress_y, stress_z;
	std::vector<double> sheer_stress_xy, sheer_stress_yz, sheer_stress_zx;
	std::vector<double> strain_x, strain_y, strain_z;
	std::vector<double> sheer_strain_xy, sheer_strain_yz, sheer_strain_zx;
	// boundary conditions
	std::vector<bool> is_boundary;
	std::vector<double> boundary_shape_function;
	std::vector<bool> is_dirichlet_vx, is_dirichlet_vy, is_dirichlet_vz, is_dirichlet_pressure;
	std::vector<double> dirichlet_vx, dirichlet_vy, dirichlet_vz, dirichlet_pressure;
	std::vector<std::vector<double>> neumann_vx, neumann_vy, neumann_vz, neumann_pressure;
	// basis function
	std::vector< std::vector<double> > element_func;
	// coloring
	std::vector<std::vector<int>> colored_elem_id;
};

struct AdjacencyMatrix{
	std::vector< std::vector<int> > idx;
	std::vector< std::vector<double> > strain_disp, stress_strain;
	std::vector< std::vector<double> > mass, stiff, damping;
};

extern SimulationParameter sim_prm;
extern StructureParameter structure;
extern AdjacencyMatrix adj_matrix;