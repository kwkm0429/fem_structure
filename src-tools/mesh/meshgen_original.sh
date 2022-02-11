#!/bin/bash

if [ $# -ne 4 ]; then
  echo "Usage: ./meshgen_origin.sh <length x> <length y> <division x> <division y>"
  exit 1;
fi

len_x=$1
len_y=$2
div_x=$3
div_y=$4
output_dir="../../data-input"
node_file="node.dat"
elem_file="elem.dat"
bc_file="bc.dat"

./bin/2d_quad_meshgen $len_x $len_y $div_x $div_y $output_dir $node_file $elem_file $bc_file