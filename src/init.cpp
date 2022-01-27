/*
* @Author: Kosuke Kawakami
* @Date:   2019-08-06 19:54:14
* @Last Modified by:   Kosuke
* @Last Modified time: 2020-08-23 11:59:19
*/
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include "parameter.h"
#include "init.h"

/*  parameter setting */
SimulationParameter sim_prm;
StructureParameter structure;
AdjacencyMatrix adj_matrix;

/**
 * @brief      Reads an input data file.
 *
 * @return     void
 */
void readNodeDataFile(){
    int loop=0, idx;
    
    std::ifstream ifsi(sim_prm.input_data_dirname+sim_prm.node_filename);   
    std::string str;
    if (ifsi.fail()) {
        std::cerr << "Failed to open " << sim_prm.node_filename << std::endl;
        exit(1);
    }
    while (std::getline(ifsi, str)) {
        std::stringstream ss(str);
        std::string item;
        std::vector<std::string> list;
        while(std::getline(ss,item,' ')&& !item.empty()){
            list.push_back(item);
        }
        if(loop==0){
            structure.num_nodes=std::stoi(list[0]);
            /* initialize fluid arryas */
            structure.x               = std::vector<double>(structure.num_nodes,0);
            structure.y               = std::vector<double>(structure.num_nodes,0);
            structure.z               = std::vector<double>(structure.num_nodes,0);
            structure.disp_x          = std::vector<double>(structure.num_nodes,0);
            structure.disp_y          = std::vector<double>(structure.num_nodes,0);
            structure.disp_z          = std::vector<double>(structure.num_nodes,0);
            structure.stress_x        = std::vector<double>(structure.num_nodes,0);
            structure.stress_y        = std::vector<double>(structure.num_nodes,0);
            structure.stress_z        = std::vector<double>(structure.num_nodes,0);
            structure.sheer_stress_xy = std::vector<double>(structure.num_nodes,0);
            structure.sheer_stress_yz = std::vector<double>(structure.num_nodes,0);
            structure.sheer_stress_zx = std::vector<double>(structure.num_nodes,0);
            structure.strain_x        = std::vector<double>(structure.num_nodes,0);
            structure.strain_y        = std::vector<double>(structure.num_nodes,0);
            structure.strain_z        = std::vector<double>(structure.num_nodes,0);
            structure.sheer_strain_xy = std::vector<double>(structure.num_nodes,0);
            structure.sheer_strain_yz = std::vector<double>(structure.num_nodes,0);
            structure.sheer_strain_zx = std::vector<double>(structure.num_nodes,0);

            structure.dirichlet_vx          = std::vector<double>(structure.num_nodes,0);
            structure.dirichlet_vy          = std::vector<double>(structure.num_nodes,0);
            structure.dirichlet_vz          = std::vector<double>(structure.num_nodes,0);
            structure.dirichlet_pressure    = std::vector<double>(structure.num_nodes,0);
            structure.neumann_vx            = std::vector<std::vector<double>>(structure.num_nodes, std::vector<double>(3,0));
            structure.neumann_vy            = std::vector<std::vector<double>>(structure.num_nodes, std::vector<double>(3,0));
            structure.neumann_vz            = std::vector<std::vector<double>>(structure.num_nodes, std::vector<double>(3,0));
            structure.neumann_pressure      = std::vector<std::vector<double>>(structure.num_nodes, std::vector<double>(3,0));
            structure.is_dirichlet_vx       = std::vector<bool>(structure.num_nodes,false);
            structure.is_dirichlet_vy       = std::vector<bool>(structure.num_nodes,false);
            structure.is_dirichlet_vz       = std::vector<bool>(structure.num_nodes,false);
            structure.is_dirichlet_pressure = std::vector<bool>(structure.num_nodes,false);
            structure.is_boundary           = std::vector<bool>(structure.num_nodes,false); 
            structure.boundary_shape_function = std::vector<double>(structure.num_nodes,0);
        }
        else if(loop>=1 && loop<=structure.num_nodes){
            /* coordinate of nodes */
            if(list.size()!=4){
                std::cerr << sim_prm.node_filename << " Error: numbers in one line" << std::endl;
                exit(1);
            }
            idx=std::stoi(list[0]);
            structure.x[idx] = std::stod(list[1]);
            structure.y[idx] = std::stod(list[2]);
            structure.z[idx] = std::stod(list[3]);
        }
        loop++;
    }
    std::cout<<"Succeeded in reading " << sim_prm.node_filename <<std::endl;
}

void readElemDataFile(){
    int i, loop=0, idx;
    
    std::ifstream ifsi(sim_prm.input_data_dirname+sim_prm.elem_filename);   
    std::string str;
    if (ifsi.fail()) {
        std::cerr << "Failed to open " << sim_prm.elem_filename << std::endl;
        exit(1);
    }
    while (std::getline(ifsi, str)) {
        std::stringstream ss(str);
        std::string item;
        std::vector<std::string> list;
        while(std::getline(ss,item,' ')&& !item.empty()){
            list.push_back(item);
        }
        if(loop==0){
            structure.num_elements         = std::stoi(list[0]);
            sim_prm.num_polygon_corner = std::stoi(list[1]);
            structure.element_node_table   = std::vector< std::vector<int> >(structure.num_elements,std::vector<int>(sim_prm.num_polygon_corner,0));
            structure.element_func         = std::vector< std::vector<double> >(structure.num_elements,std::vector<double>(sim_prm.num_polygon_corner,0));
        }else if(loop>=1 && loop<=structure.num_elements){
            if((int)(list.size())-1 != sim_prm.num_polygon_corner){
                std::cerr << sim_prm.elem_filename << " Error: numbers in one line" << std::endl;
                exit(1);
            }
            idx = std::stoi(list[0]);
            for(i=0;i<sim_prm.num_polygon_corner;i++){
                structure.element_node_table[idx][i] = std::stoi(list[i+1]);
            }
        }
        loop++;
    }
    std::cout<<"Succeeded in reading " << sim_prm.elem_filename <<std::endl;
}

void readBoundaryDataFile(){
    int i, loop = 0, node_id = 0;
    std::ifstream ifs(sim_prm.input_data_dirname + sim_prm.bc_filename);
    std::string str;
    if (ifs.fail()) {
        std::cerr << "Failed to open " << sim_prm.bc_filename << std::endl;
        exit(1);
    }
    while (std::getline(ifs, str)) {
        std::stringstream ss(str);
        std::string item;
        std::vector<std::string> list;
        while(std::getline(ss,item,' ')&& !item.empty()){
            list.push_back(item);
        }
        if(loop == 0){
            // number of nodes is written
        }else if(loop>=1 && loop<=structure.num_nodes){
            if(list.size()!=21){
                std::cout<<"Error: numbers in one line"<<std::endl;
                exit(1);
            }
            node_id = std::stoi(list[0]);

            structure.is_dirichlet_vx[node_id]        = std::stoi(list[1]);
            structure.dirichlet_vx[node_id]           = std::stod(list[2]);
            structure.is_dirichlet_vy[node_id]        = std::stoi(list[3]);
            structure.dirichlet_vy[node_id]           = std::stod(list[4]);
            structure.is_dirichlet_vz[node_id]        = std::stoi(list[5]);
            structure.dirichlet_vz[node_id]           = std::stod(list[6]);
            structure.is_dirichlet_pressure[node_id]  = std::stoi(list[7]);
            structure.dirichlet_pressure[node_id]     = std::stod(list[8]);

            for(i=0;i<3;i++){
                structure.neumann_vx[node_id][i]    = std::stod(list[9+i]);
            }
            for(i=0;i<3;i++){
                structure.neumann_vy[node_id][i]    = std::stod(list[12+i]);
            }
            for(i=0;i<3;i++){
                structure.neumann_vz[node_id][i]    = std::stod(list[15+i]);
            }
            for(i=0;i<3;i++){
                structure.neumann_pressure[node_id][i]    = std::stod(list[18+i]);
            }
            // check if the node is on boundary
            for(int i=1;i<21;i++){
                if(list[i]!="0"){
                    structure.is_boundary[node_id] = true;
                }
            }
        }
        loop++;
    }
    std::cout<<"Succeeded in reading " << sim_prm.bc_filename <<std::endl;
}

/**
 * @brief      Reads a parameter data file.
 *
 * @return     { description_of_the_return_value }
 */
void readParameterDataFile(){
    /* read from "simulation.prm" */
    std::ifstream ifsb(sim_prm.input_data_dirname + sim_prm.params_filename);
    std::string str;
    if (ifsb.fail()) {
        std::cerr << "Failed to open " << sim_prm.params_filename << std::endl;
        exit(1);
    }
    while (std::getline(ifsb, str)) {
        std::stringstream ss(str);
        std::string item;
        std::vector<std::string> list;
        while(std::getline(ss,item,' ')&& !item.empty()){
            list.push_back(item);
        }
        
        if(list[0] == "DIMENSION"){
            sim_prm.dim = std::stoi(list[2]);

        }else if(list[0] == "TIME_STEP"){
            sim_prm.max_time_step = std::stoi(list[2]);

        }else if(list[0] == "OUTPUT_INTERVAL"){
            sim_prm.output_interval = std::stoi(list[2]);

        }else if(list[0] == "DT"){
            sim_prm.dt = std::stod(list[2]);

        }else if(list[0] == "VISCOSITY"){
            structure.visc = std::stod(list[2]);

        }else if(list[0] == "FLUID_DENSITY"){
            structure.density = std::stod(list[2]);

        }else if(list[0] == "GRAVITY_X"){
            sim_prm.gravity_x = std::stod(list[2]);

        }else if(list[0] == "GRAVITY_Y"){
            sim_prm.gravity_y = std::stod(list[2]);
        
        }else if(list[0] == "GRAVITY_Z"){
            sim_prm.gravity_z = std::stod(list[2]);
        
        }else if(list[0] == "EQUATION_SOLVER"){
            sim_prm.eq_solver_opt = std::stoi(list[2]);

        }else{
            // std::cout<<list[0]<<std::endl;
        }
        structure.reynolds = structure.density / structure.visc;
    }

    std::cout<<"Succeeded in reading " << sim_prm.params_filename <<std::endl;
}

/**
 * @brief      Sets the initial state.
 */
void initField(){
    int i;
    for(i=0;i<structure.num_nodes;i++){
        if(structure.is_dirichlet_vx[i]){
        }else{
        }
        if(structure.is_dirichlet_vy[i]){
        }else{
        }
        if(structure.is_dirichlet_vz[i]){
        }else{
        }
        if(structure.is_dirichlet_pressure[i]){
        }else{
        }
    }
}

/**
 * @brief      Initializes the adjacency matrix.
 */
void initAdjMatrix(){
    int i, j, k, l, node_id1, node_id2;
    bool connected = false;

    // Adjecency matrix format used to assemble element stiffness matrix
    // adj_matrix.idx[i][j] := index of node adjuscent to i (j is order)
    adj_matrix.idx        = std::vector< std::vector<int> >(structure.num_nodes);
    adj_matrix.mass       = std::vector< std::vector<double> >(structure.num_nodes);
    adj_matrix.stiff      = std::vector< std::vector<double> >(structure.num_nodes);
    adj_matrix.damping    = std::vector< std::vector<double> >(structure.num_nodes);

    sim_prm.num_nonzero = 0;
    for(i=0;i<structure.num_elements;i++){
        for(j=0;j<sim_prm.num_polygon_corner;j++){
            node_id1 = structure.element_node_table[i][j];
            for(k=0;k<sim_prm.num_polygon_corner;k++){
                connected = false;
                node_id2 = structure.element_node_table[i][k];
                for(l=0;l<(int)adj_matrix.idx[node_id1].size();l++){
                    if(adj_matrix.idx[node_id1][l] == node_id2){
                        connected = true;
                        break;
                    }
                }
                if(connected)continue;
                sim_prm.num_nonzero++;
                adj_matrix.idx[node_id1].push_back(node_id2);
            }
        }
    }
}

void allocateAdjMatrix(){
    int i, size_of_column_nonzero;
    for(i=0;i<structure.num_nodes;i++){
        size_of_column_nonzero = adj_matrix.idx[i].size();
        adj_matrix.mass[i]        = std::vector<double>(size_of_column_nonzero, 0);
        adj_matrix.stiff[i]       = std::vector<double>(size_of_column_nonzero, 0);
        adj_matrix.damping[i]     = std::vector<double>(size_of_column_nonzero, 0);
    }
}