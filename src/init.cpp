#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include "init.h"
//#include "debug.h"

void readNodeDataFile(Sim& sim, Str& str){
    int loop=0, idx;
    
    std::ifstream ifsi(sim.input_data_dirname+sim.node_filename);   
    std::string string;
    if (ifsi.fail()) {
        std::cerr << "Failed to open " << sim.node_filename << std::endl;
        exit(1);
    }
    while (std::getline(ifsi, string)) {
        std::stringstream ss(string);
        std::string item;
        std::vector<std::string> list;
        while(std::getline(ss,item,' ')&& !item.empty()){
            list.push_back(item);
        }
        if(loop==0){
            str.num_nodes=std::stoi(list[0]);
            /* initialize structure vector */
            str.x               = std::vector<double>(str.num_nodes,0);
            str.y               = std::vector<double>(str.num_nodes,0);
            str.z               = std::vector<double>(str.num_nodes,0);
            str.disp_all        = std::vector<double>(str.num_nodes*2,0);
            str.disp_x          = std::vector<double>(str.num_nodes,0);
            str.disp_y          = std::vector<double>(str.num_nodes,0);
            str.disp_z          = std::vector<double>(str.num_nodes,0);

            str.dirichlet_dx          = std::vector<double>(str.num_nodes,0);
            str.dirichlet_dy          = std::vector<double>(str.num_nodes,0);
            str.dirichlet_dz          = std::vector<double>(str.num_nodes,0);
            str.is_dirichlet_dx       = std::vector<bool>(str.num_nodes,false);
            str.is_dirichlet_dy       = std::vector<bool>(str.num_nodes,false);
            str.is_dirichlet_dz       = std::vector<bool>(str.num_nodes,false);
            str.is_boundary           = std::vector<bool>(str.num_nodes,false); 
            str.boundary_shape_function = std::vector<double>(str.num_nodes,0);
            str.force_x               = std::vector<double>(str.num_nodes,0);
            str.force_y               = std::vector<double>(str.num_nodes,0);
            str.force_z               = std::vector<double>(str.num_nodes,0);
            str.youngs_modulus_nodes  = std::vector<double>(str.num_nodes,0);
            str.sensitivity           = std::vector<double>(str.num_nodes,0);
            str.buckling_coeff        = std::vector<double>(str.num_nodes,0);
            str.buckling_x            = std::vector< std::vector<double> >(sim.num_mode,std::vector<double>(str.num_nodes,0));
            str.buckling_y            = std::vector< std::vector<double> >(sim.num_mode,std::vector<double>(str.num_nodes,0));
            str.buckling_z            = std::vector< std::vector<double> >(sim.num_mode,std::vector<double>(str.num_nodes,0));
        }
        else if(loop>=1 && loop<=str.num_nodes){
            /* coordinate of nodes */
            if(list.size()!=4){
                std::cerr << sim.node_filename << " Error: numbers in one line" << std::endl;
                exit(1);
            }
            idx=std::stoi(list[0]);
            str.x[idx] = std::stod(list[1]);
            str.y[idx] = std::stod(list[2]);
            str.z[idx] = std::stod(list[3]);
        }
        loop++;
    }
#ifdef DEBUG
    debugPrintInfo("Input < "+sim.input_data_dirname+sim.node_filename);
    debugPrintInfo(__func__);
#endif
}

void readElemDataFile(Sim& sim, Str& str){
    int i, loop=0, idx;
    
    std::ifstream ifsi(sim.input_data_dirname+sim.elem_filename);   
    std::string string;
    if (ifsi.fail()) {
        std::cerr << "Failed to open " << sim.elem_filename << std::endl;
        exit(1);
    }
    while (std::getline(ifsi, string)) {
        std::stringstream ss(string);
        std::string item;
        std::vector<std::string> list;
        while(std::getline(ss,item,' ')&& !item.empty()){
            list.push_back(item);
        }
        if(loop==0){
            str.num_elements         = std::stoi(list[0]);
            sim.num_polygon_corner = std::stoi(list[1]);
            str.element_node_table   = std::vector< std::vector<int> >(str.num_elements,std::vector<int>(sim.num_polygon_corner,0));
            str.element_func         = std::vector< std::vector<double> >(str.num_elements,std::vector<double>(sim.num_polygon_corner,0));

            str.stress_x        = std::vector<double>(str.num_elements,0);
            str.stress_y        = std::vector<double>(str.num_elements,0);
            str.stress_z        = std::vector<double>(str.num_elements,0);
            str.sheer_stress_xy = std::vector<double>(str.num_elements,0);
            str.sheer_stress_yz = std::vector<double>(str.num_elements,0);
            str.sheer_stress_zx = std::vector<double>(str.num_elements,0);
            str.strain_x        = std::vector<double>(str.num_elements,0);
            str.strain_y        = std::vector<double>(str.num_elements,0);
            str.strain_z        = std::vector<double>(str.num_elements,0);
            str.sheer_strain_xy = std::vector<double>(str.num_elements,0);
            str.sheer_strain_yz = std::vector<double>(str.num_elements,0);
            str.sheer_strain_zx = std::vector<double>(str.num_elements,0);

        }else if(loop>=1 && loop<=str.num_elements){
            if((int)(list.size())-1 != sim.num_polygon_corner){
                std::cerr << sim.elem_filename << " Error: numbers in one line" << std::endl;
                exit(1);
            }
            idx = std::stoi(list[0]);
            for(i=0;i<sim.num_polygon_corner;i++){
                str.element_node_table[idx][i] = std::stoi(list[i+1]);
            }
        }
        loop++;
    }
#ifdef DEBUG
    debugPrintInfo("Input < "+sim.input_data_dirname+sim.elem_filename);
    debugPrintInfo(__func__);
#endif
}

void readBoundaryDataFile(Sim& sim, Str& str){
    int i, loop = 0, node_id = 0, list_num=0;
    std::ifstream ifs(sim.input_data_dirname + sim.bc_filename);
    std::string string;
    if (ifs.fail()) {
        std::cerr << "Failed to open " << sim.bc_filename << std::endl;
        exit(1);
    }
    while (std::getline(ifs, string)) {
        std::stringstream ss(string);
        std::string item;
        std::vector<std::string> list;
        while(std::getline(ss,item,' ')&& !item.empty()){
            list.push_back(item);
        }
        if(loop == 0){
            // number of nodes is written
        }else if(loop>=1 && loop<=str.num_nodes){
            if(list.size()!=11){
                std::cout<<"Error: numbers in one line"<<std::endl;
                exit(1);
            }
            node_id = std::stoi(list[0]);
            list_num = std::stoi(list[1]);

            str.is_dirichlet_dx[node_id]        = std::stoi(list[2]);
            str.dirichlet_dx[node_id]           = std::stod(list[3]);
            str.is_dirichlet_dy[node_id]        = std::stoi(list[4]);
            str.dirichlet_dy[node_id]           = std::stod(list[5]);
            str.is_dirichlet_dz[node_id]        = std::stoi(list[6]);
            str.dirichlet_dz[node_id]           = std::stod(list[7]);
            str.force_x[node_id] = std::stod(list[8]);
            str.force_y[node_id] = std::stod(list[9]);
            str.force_z[node_id] = std::stod(list[10]);

            // check if the node is on boundary
            for(int i=1;i<11;i++){
                if(list[i]!="0"){
                    str.is_boundary[node_id] = true;
                }
            }
        }
        loop++;
    }
#ifdef DEBUG
    debugPrintInfo("Input < "+sim.input_data_dirname + sim.bc_filename);
    debugPrintInfo(__func__);
#endif
}

void readParameterDataFile(Sim& sim, Str& str){
    /* read from "simulation.prm" */
    std::ifstream ifsb(sim.input_data_dirname + sim.params_filename);
    std::string string;
    if (ifsb.fail()) {
        std::cerr << "Failed to open " << sim.params_filename << std::endl;
        exit(1);
    }
    while (std::getline(ifsb, string)) {
        std::stringstream ss(string);
        std::string item;
        std::vector<std::string> list;
        while(std::getline(ss,item,' ')&& !item.empty()){
            list.push_back(item);
        }

        if(list[0] == "YOUNGS_MODULUS"){
            str.youngs_modulus = std::stod(list[2]);
        
        }else if(list[0] == "POISSON_RATIO"){
            str.poisson_ratio = std::stod(list[2]);

        }else if(list[0] == "DIMENSION"){
            sim.dim = std::stoi(list[2]);

        }else if(list[0] == "TIME_STEP"){
            sim.max_time_step = std::stoi(list[2]);

        }else if(list[0] == "OUTPUT_INTERVAL"){
            sim.output_interval = std::stoi(list[2]);

        }else if(list[0] == "DT"){
            sim.dt = std::stod(list[2]);

        }else if(list[0] == "DENSITY"){
            str.density = std::stod(list[2]);

        }else if(list[0] == "GRAVITY_X"){
            sim.gravity_x = std::stod(list[2]);

        }else if(list[0] == "GRAVITY_Y"){
            sim.gravity_y = std::stod(list[2]);
        
        }else if(list[0] == "GRAVITY_Z"){
            sim.gravity_z = std::stod(list[2]);
        
        }else if(list[0] == "EQUATION_SOLVER"){
            sim.eq_solver_opt = std::stoi(list[2]);

        }else if(list[0] == "THICKNESS"){
            str.thickness = std::stod(list[2]);
        
        }else if(list[0] == "ID"){
            sim.id = std::stoi(list[2]);
        
        }else if(list[0] == "NUM_MODE"){
            sim.num_mode = std::stoi(list[2]);
        
        }else{
            // std::cout<<list[0]<<std::endl;
        }
    }
#ifdef DEBUG
    debugPrintInfo("Input < "+sim.input_data_dirname + sim.params_filename);
    debugPrintInfo(__func__);
#endif
}

void initField(Sim& sim, Str& str){
    int i;
    for(i=0;i<str.num_nodes;i++){
        if(str.is_dirichlet_dx[i]){
        }else{
        }
        if(str.is_dirichlet_dy[i]){
        }else{
        }
        if(str.is_dirichlet_dz[i]){
        }else{
        }
    }
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void initAdjMatrix(Sim& sim, Str& str, AdjMatrix& adj_mat){
    int i, j, k, l, node_id1, node_id2;
    bool connected = false;

    // Adjecency matrix format used to assemble element stiffness matrix
    // adj_mat.idx[i][j] := index of node adjuscent to i (j is order)
    adj_mat.idx        = std::vector< std::vector<int> >(str.num_nodes);
    adj_mat.stiff      = std::vector< std::vector<double> >(str.num_nodes*sim.dim);
    adj_mat.stiff_geo  = std::vector< std::vector<double> >(str.num_nodes*sim.dim);

    sim.num_nonzero = 0;
    for(i=0;i<str.num_elements;i++){
        for(j=0;j<sim.num_polygon_corner;j++){
            node_id1 = str.element_node_table[i][j];
            for(k=0;k<sim.num_polygon_corner;k++){
                connected = false;
                node_id2 = str.element_node_table[i][k];
                for(l=0;l<(int)adj_mat.idx[node_id1].size();l++){
                    if(adj_mat.idx[node_id1][l] == node_id2){
                        connected = true;
                        break;
                    }
                }
                if(connected)continue;
                sim.num_nonzero++;
                adj_mat.idx[node_id1].push_back(node_id2);
            }
        }
    }
    sim.num_nonzero*=4;
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void allocateAdjMatrix(Sim& sim, Str& str, AdjMatrix& adj_mat){
    int i, size_of_column_nonzero;
    for(i=0;i<str.num_nodes;i++){
        size_of_column_nonzero = adj_mat.idx[i].size()*sim.dim;
        adj_mat.stiff[i*2]      = std::vector<double>(size_of_column_nonzero, 0);
        adj_mat.stiff[i*2+1]    = std::vector<double>(size_of_column_nonzero, 0);
        adj_mat.stiff_geo[i*2]   = std::vector<double>(size_of_column_nonzero, 0);
        adj_mat.stiff_geo[i*2+1] = std::vector<double>(size_of_column_nonzero, 0);
    }
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}