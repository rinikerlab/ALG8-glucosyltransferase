c=$1
let s=c+1
input=ALG8_pretransfer_H40d_DSG_prod$c
output=ALG8_pretransfer_H40d_DSG_prod$s
pmemd.cuda -O -i ./mdh.in -p ./ALG8_pretransfer_H40d_DSG.prmtop -c $input.rst -o $output.out -r $output.rst -ref $input.rst -x $output.nc


