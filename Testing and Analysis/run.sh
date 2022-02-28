#!/usr/bin/bash

echo "Starting simulation"
./../build/MCScattering
wait
echo "Starting ROI analysis"
path=$(<../build/outputpath.txt)
path="${path}/Blurred Images/"
python /mnt/c/Users/adam/Documents/Code/MontePython/Testing\ and\ Analysis/raw.py -i "$path"