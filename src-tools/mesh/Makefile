CXX     := g++  

all: 2d_quad_meshgen 3d_msh2dat 2d_msh2dat write_geo 3d_msh2dat_str clean

clean:
	rm ./*.o

# make 2d mesh generator
2d_quad_meshgen: 2d_quad_meshgen.o
	$(CXX) -o bin/2d_quad_meshgen 2d_quad_meshgen.o

2d_quad_meshgen.o:
	$(CXX) -std=c++11 -o 2d_quad_meshgen.o -c src/2d_quad_meshgen.cpp

# make 3d mesh converter
3d_msh2dat: 3d_msh2dat.o
	$(CXX) -o bin/3d_msh2dat 3d_msh2dat.o

3d_msh2dat.o:
	$(CXX) -std=c++11 -o 3d_msh2dat.o -c src/3d_msh2dat.cpp

# make 2d mesh converter
2d_msh2dat: 2d_msh2dat.o
	$(CXX) -o bin/2d_msh2dat 2d_msh2dat.o

2d_msh2dat.o:
	$(CXX) -std=c++11 -o 2d_msh2dat.o -c src/2d_msh2dat.cpp

write_geo: write_geo.o
	$(CXX) -o bin/write_geo write_geo.o

write_geo.o:
	$(CXX) -std=c++11 -o write_geo.o -c src/write_geo.cpp

3d_msh2dat_str: 3d_msh2dat_str.o
	$(CXX) -o bin/3d_msh2dat_str 3d_msh2dat_str.o

3d_msh2dat_str.o:
	$(CXX) -std=c++11 -o 3d_msh2dat_str.o -c src/3d_msh2dat_str.cpp
