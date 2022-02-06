#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

int dim = 2;
int div_num = 20;
std::string geo_filename;

void outputGeoFile(){
	std::ofstream ofs(geo_filename);
    if(dim == 2){
        ofs<<"\
Point(1) = {0, 0, 0, 1.0};\n\
Point(2) = {1, 0, 0, 1.0};\n\
Point(3) = {1, 1, 0, 1.0};\n\
Point(4) = {0, 1, 0, 1.0};\n\
Line(1) = {1, 2};\n\
Line(2) = {2, 3};\n\
Line(3) = {3, 4};\n\
Line(4) = {4, 1};\n\
Line Loop(1) = {4, 1, 2, 3};\n\
Plane Surface(1) = {1};\n\
Transfinite Line {2, 3, 4, 1} = " + std::to_string(div_num+1) + " Using Progression 1;\n\
Transfinite Surface {1};\n\
Recombine Surface {1};\n\
Physical Line(\"movingWall\") = {3};\n\
Physical Line(\"wall\") = {4, 1, 2};\n\
Physical Surface(\"fluid\") = {1};\n\
        ";
    }else if(dim == 3){
        ofs<<"\
Point(1) = {0, 0, 0, 1.0};\n\
Point(2) = {1, 0, 0, 1.0};\n\
Point(3) = {1, 1, 0, 1.0};\n\
Point(4) = {0, 1, 0, 1.0};\n\
Point(5) = {0, 0, 1, 1.0};\n\
Point(6) = {1, 0, 1, 1.0};\n\
Point(7) = {1, 1, 1, 1.0};\n\
Point(8) = {0, 1, 1, 1.0};\n\
Point(9) = {1, 1, 1, 1.0};\n\
Line(1) = {1, 2};\n\
Line(2) = {2, 3};\n\
Line(3) = {3, 4};\n\
Line(4) = {4, 1};\n\
Line(5) = {1, 5};\n\
Line(6) = {2, 6};\n\
Line(7) = {3, 7};\n\
Line(8) = {4, 8};\n\
Line(9) = {8, 7};\n\
Line(10) = {7, 6};\n\
Line(11) = {6, 5};\n\
Line(12) = {5, 8};\n\
Line Loop(1) = {1, 2, 3, 4};\n\
Surface(1) = {1};\n\
Line Loop(3) = {3, 8, 9, -7};\n\
Surface(2) = {3};\n\
Line Loop(5) = {4, 5, 12, -8};\n\
Surface(3) = {5};\n\
Line Loop(7) = {7, 10, -6, 2};\n\
Surface(4) = {7};\n\
Line Loop(9) = {1, 6, 11, -5};\n\
Surface(5) = {9};\n\
Line Loop(11) = {9, 10, 11, 12};\n\
Surface(6) = {11};\n\
Surface Loop(1) = {6, 2, 4, 5, 1, 3};\n\
Volume(1) = {1};\n\
Transfinite Line \"*\" = " + std::to_string(div_num+1) + " Using Bump 1.0;\n\
Transfinite Surface \"*\";\n\
Transfinite Volume \"*\";\n\
Physical Surface(\"movingWall\") = {6};\n\
Physical Surface(\"wall\") = {1, 2, 3, 4, 5};\n\
Physical Volume(\"fluid\") = {1};\n\
        ";
    }else{
        std::cout<<"Dimension Error: <dimension> has to be 2 or 3"<<std::endl;
        exit(1);
    }

	std::cout<<"Finish outputGeoFile()"<<std::endl;
}

int main(int argc, char* argv[]){
    if(argc != 4){
        fprintf(stderr,"Usage: ./write_geo <geo file> <dimension> <division number>\n");
        exit(1);
    }

    geo_filename = argv[1];
    dim = atoi(argv[2]);
    div_num = atoi(argv[3]);

	outputGeoFile();
}
