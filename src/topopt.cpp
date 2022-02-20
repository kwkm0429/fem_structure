#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

#include "topopt.h"
#include "structure_solver.h"
#include "eigen_solver.h"
#include "output.h"
#include "init.h"
//#include "debug.h"
#include "fem_matrix.h"

void readTopOptDataFile(TopOpt& top, Sim& sim){
    /* read from "simulation.prm" */
    std::ifstream ifsb(sim.input_data_dirname + top.params_filename);
    std::string string;
    if (ifsb.fail()) {
        std::cerr << "Failed to open " << top.params_filename << std::endl;
        exit(1);
    }
    while (std::getline(ifsb, string)) {
        std::stringstream ss(string);
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
    debugPrintInfo("Input < "+sim.input_data_dirname + top.params_filename);
    debugPrintInfo(__func__);
#endif
}

void initTopOpt(TopOpt& top, Str& str){
	int i;
    top.rho = std::vector<double>(str.num_nodes,0);
    top.rho_new = std::vector<double>(str.num_nodes,0);
    top.sens = std::vector<double>(str.num_nodes,0);
	for(i=0;i<str.num_nodes;i++){
		top.rho[i] = top.rho_init;
	}
    top.comp = top.comp_prev = 0;
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void calcVolume(TopOpt& top, Str& str){
    int i, j, node_id;
    double vol_domain = 0;
    top.vol_sum = 0;
    top.vol_frac = 0;
    for(i=0;i<str.num_elements;i++){
        for(j=0;j<4;j++){
            node_id = str.element_node_table[i][j];
            top.vol_sum += top.rho[node_id] * str.element_func[i][j];
            vol_domain += str.element_func[i][j];
        }
    }
    top.vol_frac = top.vol_sum / vol_domain;
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void optOCmethod(TopOpt& top){
    int i;
    int rsize=(int)(top.rho.size());
    double lambda1=0, lambda2=1e4, lmid=0;
    double mvlmt=0.15;
    double eta=0.5;
    double mean;
    std::vector<double> M = std::vector<double>(rsize,0);
    while((lambda2-lambda1)/(lambda1+lambda2)>1e-3){
        lmid = (lambda2+lambda1)*0.5;
        for(i=0;i<rsize;i++){
            if(top.sens[i]==0)M[i]=0;
            else M[i] = top.rho[i] * std::pow((-top.sens[i]/lmid), eta);
            top.rho_new[i] = std::max(std::max(std::min(std::min(M[i], top.rho[i]+mvlmt), 1.0), top.rho[i]-mvlmt), 1e-9);
        }
        mean = 0;
        for(i=0;i<rsize;i++){
            mean += top.rho_new[i];
        }
        mean /= rsize;
        if(mean-top.vol_max>0){
            lambda1 = lmid;
        }else{
            lambda2 = lmid;
        }
    }
    for(i=0;i<(int)(top.rho.size());i++){
        top.rho[i] = top.rho_new[i];
    }
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void updateYoungsModulus(TopOpt& top, Str& str){
    int i;
    for(i=0;i<str.num_nodes;i++){
        str.youngs_modulus_nodes[i] = (top.E0-top.Emin)*std::pow(top.rho[i], top.pow)+top.Emin;
    }
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}

void exeTopOpt(TopOpt& top, Sim& sim, Str& str, AdjMatrix& adj_mat){
    int i;
    std::cout<<"---------- Topology Optimization Start ----------"<<std::endl;
    // init topopt
    readTopOptDataFile(top, sim);
    initTopOpt(top, str);
    outputDensityVtkFile(0, top, sim, str);
    outputTopOptDataFile(top, sim);
    
    for(i=0;i<top.itr_max;i++){
        updateYoungsModulus(top, str);
        exeStaticAnalysis(sim, str, adj_mat);
        calcVolume(top, str);
        calcCompliance(top.comp);
        if(std::abs((top.comp - top.comp_prev)*(top.comp - top.comp_prev))<top.comp_conv){// convergence check
            std::cout<<"converged"<<std::endl;
            break;
        }
        top.comp_prev = top.comp;
        calcElementMatrix2Dquad(sim, str, adj_mat);
        calcSensitivity(top, str);
        optOCmethod(top);
        outputDensityVtkFile(i+1, top, sim, str);
        std::cout.setf(std::ios::left, std::ios::adjustfield);
        std::cout<<" Step: "<<std::setw(4)<<i+1;
        std::cout<<" Volume: "<<std::setw(10)<<top.vol_frac;
        std::cout<<" Compliance: "<<std::setw(10)<<top.comp<<std::endl;
    }
#ifdef DEBUG
    debugPrintInfo(__func__);
#endif
}