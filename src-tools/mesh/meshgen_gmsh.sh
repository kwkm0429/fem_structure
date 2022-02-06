#!/bin/bash

if [ $# -ne 4 ]; then
  echo "Usage: ./meshgen_gmsh.sh <dimension> <division> <.geo file> <.msh file>"
  exit 1;
fi

dim=$1
div=$2
geo_file=$3
msh_file=$4
output_dir="../../data-input"
node_file="node.dat"
elem_file="elem.dat"
bc_file="bc.dat"

./bin/write_geo $geo_file $dim $div
gmsh -3 $geo_file -format msh2 -o $msh_file
if [ $dim -eq 2 ]; then
	./bin/2d_msh2dat $msh_file $output_dir $node_file $elem_file $bc_file
elif [ $dim -eq 3 ]; then
	./bin/3d_msh2dat $msh_file $output_dir $node_file $elem_file $bc_file
else 
	echo "input <dimension> is "$dim
	echo "Warning: <dimension> has to be 2 or 3"
fi