#!/bin/bash

if [ $# -ne 2 ]; then
  echo "Usage: ./gmsh_geo2msh.sh <geo_file> <msh_file>"
  exit 1;
fi

geo_file=$1
msh_file=$2

gmsh -3 $geo_file -format msh2 -o $msh_file