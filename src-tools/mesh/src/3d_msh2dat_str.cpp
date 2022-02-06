/*
* @Author: Kosuke
* @Date:   2020-06-26 13:27:46
* @Last Modified by:   Kosuke
* @Last Modified time: 2021-08-28 23:34:47
*/
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

std::string mesh_filename  = "";
std::string output_dirname = "";
std::string node_filename  = "";
std::string elem_filename  = "";

int num_nodes, num_elements;
std::vector<double> node_x, node_y, node_z;
std::vector<std::vector<int>> element_node_table;

struct physics{
	int dim;
	std::string name;
};

std::vector<physics> phy = std::vector<physics>();

bool is_equal(double a, double b){
    if(std::abs(a-b)<1e-8)return true;
    else return false;
}

void setVectorOfNode(int num){
	node_x = std::vector<double>(num,0);
	node_y = std::vector<double>(num,0);
	node_z = std::vector<double>(num,0);
}

void readMeshDataFile(){
	bool is_physics_name = false;
	bool is_nodes = false;
	bool is_elements = false;
	int loop = 0;

	std::ifstream ifsi(mesh_filename);
    std::string str;
	if (ifsi.fail()) {
        std::cerr << "Failed to open" << mesh_filename << std::endl;
        exit(1);
    }
    while (std::getline(ifsi, str)) {
        std::stringstream ss(str);
        std::string item;
        std::vector<std::string> list;
        while(std::getline(ss,item,' ')&& !item.empty()){
            list.push_back(item);
        }
        if(list.size()==1)std::cout<<list[0]<<std::endl;
        // PhysicalName
        if(list[0].size() >= 13 && list[0].substr(1, 13) == "PhysicalNames"){
        	is_physics_name = true;
        	loop = 0;
        	continue;
        }
        if(is_physics_name){
        	if(loop == 0){
        		int num = std::stoi(list[0]);
        		phy = std::vector<physics>(num);
        	}else if(loop>0 && loop<=(int)(phy.size())){
        		int id = loop - 1;
        		phy[id].dim = std::stoi(list[0]);
                int len = list[2].size();
                phy[id].name = list[2].substr(1, len-3);
        	}
        	loop++;
        }
        if(list[0].size() >= 16 && list[0].substr(1, 16) == "EndPhysicalNames"){
        	is_physics_name = false;
        	loop = 0;
        }
        // Node
        if(list[0].size() >= 5 && list[0].substr(1, 5) == "Nodes"){
        	is_nodes = true;
        	loop = 0;
        	continue;
        }
        if(is_nodes){
        	if(loop == 0){
        		num_nodes = std::stoi(list[0]);
        		setVectorOfNode(num_nodes);
        	}else if(loop>0 && loop<=num_nodes){
        		int id = loop - 1;
        		node_x[id] = std::stod(list[1]);
        		node_y[id] = std::stod(list[2]);
        		node_z[id] = std::stod(list[3]);
        	}
        	loop++;
        }
        if(list[0].size() >= 8 && list[0].substr(1, 8) == "EndNodes"){
        	is_nodes = false;
        	loop = 0;
        }
        // Elements
        if(list[0].size() >= 8 && list[0].substr(1, 8) == "Elements"){
        	is_elements = true;
        	loop = 0;
        	continue;
        }
        if(is_elements){
        	if(loop == 0){
        		num_elements = std::stoi(list[0]);
        	}else if(loop>0 && loop<= num_elements){
        		int num = std::stoi(list[1]);
        		if(num == 4 || num == 3){
        			// element
        			int id1 = std::stoi(list[5]) - 1;
        			int id2 = std::stoi(list[6]) - 1;
        			int id3 = std::stoi(list[7]) - 1;
        			int id4 = std::stoi(list[8]) - 1;
        			element_node_table.push_back(std::vector<int>{id1, id2, id3, id4});
        		}
        	}
            loop++;
        }
        if(list[0].size() >= 11 && list[0].substr(1, 11) == "EndElements"){
        	is_elements = false;
        	loop = 0;
        }
    }
    std::cout<<"Succeeded readNodeDataFile()"<<std::endl;
}

void outputNodeFile(){
	int i;
	std::ofstream ofs(output_dirname + "/" + node_filename);
	ofs<<num_nodes<<std::endl;
	for(i=0;i<num_nodes;i++){
		ofs<<node_x[i]<<" "<<node_y[i]<<" "<<node_z[i]<<std::endl;
	}
	std::cout<<"Finish outputNodeFile()"<<std::endl;
    std::cout<<"Output > "<<output_dirname+"/"+node_filename<<std::endl;
}

void outputElemFile(){
	int i,j;
    int size = element_node_table.size();
	std::ofstream ofs(output_dirname + "/" + elem_filename);
	ofs<<size<<" "<<4<<std::endl;
	for(i=0;i<size;i++){
		for(j=0;j<4;j++){
			ofs<<element_node_table[i][j]<<" ";
		}
		ofs<<std::endl;
	}
	std::cout<<"Finish outputElemFile()"<<std::endl;
    std::cout<<"Output > "<<output_dirname+"/"+elem_filename<<std::endl;
}

int main(int argc, char* argv[]){
    if(argc != 5){
        fprintf(stderr,"Usage: ./3d_msh2dat_str <msh file> <output directory> <node filename> <elem file name>\n");
        exit(1);
    }

    mesh_filename = argv[1];
    output_dirname = argv[2];
    node_filename = argv[3];
    elem_filename = argv[4];

	readMeshDataFile();
	outputNodeFile();
	outputElemFile();
}