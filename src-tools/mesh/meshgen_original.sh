#!/bin/bash

if [ $# -ne 0 ]; then
  echo "Usage: ./meshgen_origin.sh"
  exit 1;
fi

output_dir="../../data-input"
node_file="node.dat"
elem_file="elem.dat"
bc_file="bc.dat"

./bin/2d_quad_meshgen $output_dir $node_file $elem_file $bc_file