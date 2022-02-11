#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <stdlib.h>

// default parameters for structure mesh
double LENGTH_X = 2.0;
double LENGTH_Y = 1.0;
double DX = 0.2;
double DY = 0.2;
int SIZE_X=LENGTH_X*(1/DX) + 1;
int SIZE_Y=LENGTH_Y*(1/DY) + 1;
int NUMBER_OF_NODES=SIZE_Y * SIZE_X;
int NUMBER_OF_ELEMENTS=(SIZE_Y-1) * (SIZE_X-1);
int NUM_NODES_PER_ELEM=4;

// names of input files
std::string mesh_filename  = "";
std::string output_dirname = "";
std::string node_filename  = "";
std::string elem_filename  = "";
std::string bc_filename    = "";

double pi = 3.141592653589793;

std::vector<double> NodeX, NodeY, NodeZ;
std::vector<double> DirichletVx, DirichletVy, DirichletVz;
std::vector<std::vector<double>> NeumannVx, NeumannVy, NeumannVz;
std::vector<bool> isDirichletVx, isDirichletVy, isDirichletVz;
std::vector<double> ForceX, ForceY, ForceZ;
std::vector< std::vector<int> > ElementNodeTable;

void setNodePosition();
void setElementNodeTable();
void setBoundaryCondition();
bool equal(double, double);
void initFluidVector();
void makeNodeFile();
void makeElemFile();
void makeBoundaryFile();

int main(int argc, char* argv[]){

	if(argc != 9){
        fprintf(stderr,"Usage: ./2d_quad_meshgen <length x> <length y> <division x> <division y> <output directory> <node filename> <elem file name> <bc file name>\n");
        exit(1);
    }

    LENGTH_X = atof(argv[1]);
	LENGTH_Y = atof(argv[2]);
	SIZE_X = atof(argv[3])+1;
	SIZE_Y = atof(argv[4])+1;
	DX = LENGTH_X / (SIZE_X - 1);
	DY = LENGTH_Y / (SIZE_Y - 1);
	NUMBER_OF_NODES=SIZE_Y * SIZE_X;
	NUMBER_OF_ELEMENTS=(SIZE_Y-1) * (SIZE_X-1);
	NUM_NODES_PER_ELEM=4;

    output_dirname = argv[5];
    node_filename = argv[6];
    elem_filename = argv[7];
    bc_filename = argv[8];

	initFluidVector();
	setNodePosition();
	setElementNodeTable();
	setBoundaryCondition();
	makeNodeFile();
	makeElemFile();
	makeBoundaryFile();
}

void setNodePosition(){
	int i, j, nodeIdx;
	for(i=0;i<SIZE_Y;i++){
		for(j=0;j<SIZE_X;j++){
			nodeIdx=i*SIZE_X+j;
			/* set points */
			NodeX[nodeIdx]=j*DX;
			NodeY[nodeIdx]=i*DY;
			NodeZ[nodeIdx]=0;
		}
	}
}

void setElementNodeTable(){
	int i, ei, ej;
	for(i=0;i<NUMBER_OF_ELEMENTS;i++){
		ei=i/(SIZE_X-1);
		ej=i%(SIZE_X-1);
		
		ElementNodeTable[i][0]=ei*SIZE_X+ej;
		ElementNodeTable[i][1]=ei*SIZE_X+ej+1;
		ElementNodeTable[i][2]=(ei+1)*SIZE_X+ej+1;
		ElementNodeTable[i][3]=(ei+1)*SIZE_X+ej;
	}
}

void setBoundaryCondition(){
	int i;
	for(i=0;i<NUMBER_OF_NODES;i++){
		//* Cantilever
		if(equal(NodeX[i],0)){
			isDirichletVx[i]=true; DirichletVx[i]=0;
			isDirichletVy[i]=true; DirichletVy[i]=0;
			isDirichletVz[i]=true; DirichletVz[i]=0;
		}
		if(equal(NodeX[i], LENGTH_X)){
			ForceX[i]=0;
			ForceY[i]=-10000;
			ForceZ[i]=0;
		}
		//*/
	}
}

bool equal(double a, double b){
	if(std::abs(a-b)<1e-9)return true;
	else return false;
}

void initFluidVector(){
	NodeX=std::vector<double>(NUMBER_OF_NODES,0);
	NodeY=std::vector<double>(NUMBER_OF_NODES,0);
	NodeZ=std::vector<double>(NUMBER_OF_NODES,0);
	ElementNodeTable=std::vector< std::vector<int> >(NUMBER_OF_ELEMENTS,std::vector<int>(NUM_NODES_PER_ELEM,0));

	DirichletVx=std::vector<double>(NUMBER_OF_NODES,0);
	DirichletVy=std::vector<double>(NUMBER_OF_NODES,0);
	DirichletVz=std::vector<double>(NUMBER_OF_NODES,0);
	isDirichletVx=std::vector<bool>(NUMBER_OF_NODES,0);
	isDirichletVy=std::vector<bool>(NUMBER_OF_NODES,0);
	isDirichletVz=std::vector<bool>(NUMBER_OF_NODES,0);
	ForceX=std::vector<double>(NUMBER_OF_NODES,0);
	ForceY=std::vector<double>(NUMBER_OF_NODES,0);
	ForceZ=std::vector<double>(NUMBER_OF_NODES,0);
}

void makeNodeFile(){
	int i;
	std::ofstream ofs(output_dirname + "/" + node_filename);
	ofs<<NUMBER_OF_NODES<<std::endl;
	for(i=0;i<NUMBER_OF_NODES;i++){
		ofs<<i<<" "<<NodeX[i]<<" "<<NodeY[i]<<" "<<NodeZ[i]<<std::endl;
	}
	ofs.close();
	std::cout<<"Finish makeNodeFile()"<<std::endl;
	std::cout<<"Output > "<<output_dirname + "/" + node_filename<<std::endl;
}

void makeElemFile(){
	int i,j;
	std::ofstream ofs(output_dirname + "/" + elem_filename);

	ofs<<NUMBER_OF_ELEMENTS<<" "<<NUM_NODES_PER_ELEM<<std::endl;
	for(i=0;i<NUMBER_OF_ELEMENTS;i++){
		ofs<<i<<" ";
		for(j=0;j<NUM_NODES_PER_ELEM;j++){
			ofs<<ElementNodeTable[i][j];
			if(j!=NUM_NODES_PER_ELEM-1)ofs<<" ";
			else ofs<<std::endl;
		}
	}
	ofs.close();
	std::cout<<"Finish makeElemFile()"<<std::endl;
	std::cout<<"Output > "<<output_dirname + "/" + elem_filename<<std::endl;
}

void makeBoundaryFile(){
	int i, j;
	std::ofstream ofs(output_dirname + "/" + bc_filename);
	ofs<<NUMBER_OF_NODES<<" "<<10<<std::endl;
	for(i=0;i<NUMBER_OF_NODES;i++){
		ofs<<i<<" "<<10;
		ofs<<" "<<isDirichletVx[i]<<" "<<DirichletVx[i]<<
		     " "<<isDirichletVy[i]<<" "<<DirichletVy[i]<<
		     " "<<isDirichletVz[i]<<" "<<DirichletVz[i]<<
		     " "<<ForceX[i]<<" "<<ForceY[i]<<" "<<ForceZ[i]<<std::endl;
	}
	ofs.close();
	std::cout<<"Finish makeBoundaryFile()"<<std::endl;
	std::cout<<"Output > "<<output_dirname + "/" + bc_filename<<std::endl;
}