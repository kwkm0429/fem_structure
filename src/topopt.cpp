#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include "topopt.h"
#include "parameter.h"
#include "structure_solver.h"
#include "eigen_solver.h"
#include "output.h"
#include "init.h"
#include "debug.h"

TopOptParameter top;

void readTopOptDataFile(){
    /* read from "simulation.prm" */
    std::ifstream ifsb(sim_prm.input_data_dirname + top.params_filename);
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

        if(list[0] == "VOL_MAX"){
            top.vol_max = std::stod(list[2]);

        }else if(list[0] == "RHO_INIT"){
            top.rho_init = std::stod(list[2]);
        
        }else if(list[0] == "E0"){
            top.E0 = std::stod(list[2]);

        }else if(list[0] == "Emin"){
            top.Emin = std::stod(list[2]);

        }else if(list[0] == "POW"){
            top.pow = std::stoi(list[2]);

        }else if(list[0] == "ITR_MAX"){
            top.itr_max = std::stoi(list[2]);

        }else{
            // std::cout<<list[0]<<std::endl;
        }
    }
#ifdef DEBUG
    debugPrintInfo("Input < "+sim_prm.input_data_dirname + top.params_filename);
    debugPrintInfo(__func__);
#endif
}

void initTopOpt(){
	int i;
    top.rho = std::vector<double>(structure.num_nodes,0);
    top.sens = std::vector<double>(structure.num_nodes,0);
	for(i=0;i<structure.num_nodes;i++){
		top.rho[i] = top.rho_init;
	}
    top.comp = top.comp_prev = 0;
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void calcVolume(){
    int i, j, node_id;
    double vol_domain = 0;
    top.vol_sum = 0;
    top.vol_frac = 0;
    for(i=0;i<structure.num_elements;i++){
        for(j=0;j<4;j++){
            node_id = structure.element_node_table[i][j];
            top.vol_sum += top.rho[node_id] * structure.element_func[i][j];
            vol_domain += structure.element_func[i][j];
        }
    }
    top.vol_frac = top.vol_sum / vol_domain;
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void optOCmethod(){
    int i, rsize=top.rho.size();
    double lambda1=0, lambda2=1e4, lmid=0;
    double mvlmt = 0.15;
    double eta=0.5;
    double mean=0;
    std::vector<double> M = std::vector<double>(rsize,0);
    while((lambda2-lambda1)/(lambda1+lambda2)>1e-3){
        lmid=(lambda2+lambda1)*0.5;
        for(i=0;i<rsize;i++){
            M[i] = top.rho[i] * std::pow((-top.sens[i]/lmid), eta);
            top.rho[i] = std::max(std::max(std::min(std::min(M[i], top.rho[i]+mvlmt), 1.0), top.rho[i]-mvlmt), 0.0);
        }
        for(i=0;i<rsize;i++){
            mean += top.rho[i];
        }
        mean /= rsize;
        if(mean-top.vol_max>0){
            lambda1 = lmid;
        }else{
            lambda2 = lmid;
        }
    }
}

void updateYoungsModulus(){
    int i;
    for(i=0;i<structure.num_nodes;i++){
        structure.youngs_modulus_nodes[i] = (top.E0-top.Emin)*std::pow(top.rho[i], top.pow)+top.Emin;
    }
}

void exeTopOpt(){
    int i;
    std::cout<<"---------- Topology Optimization Start ----------"<<std::endl;
    // init topopt
    readTopOptDataFile();
    initTopOpt();
    outputDensityVtkFile(0, top);
    outputTopOptDataFile(top);
    
    for(i=0;i<top.itr_max;i++){
        std::cout<<"---------- Topology Optimization Step "<<i+1<<" ----------"<<std::endl;
        updateYoungsModulus();
        exeStaticAnalysis();
        calcVolume();
        calcCompliance(top.comp);
        if(std::abs((top.comp - top.comp_prev)*(top.comp - top.comp_prev))<top.comp_conv){// convergence check
            std::cout<<"converged"<<std::endl;
            break;
        }
        top.comp_prev = top.comp;
        calcSensitivity(top);
        optOCmethod();
        outputDensityVtkFile(i+1, top);
        std::cout<<"Volume: "<<top.vol_frac<<std::endl;
        std::cout<<"Compliance: "<<top.comp<<std::endl;
    }
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}