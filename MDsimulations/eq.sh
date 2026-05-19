# equilibration
pmemd.cuda -O -i ./input_files/eq1.in -p ./ALG8_pretransfer_H40d_DSG.prmtop -c ./ALG8_pretransfer_H40d_DSG_min2.rst -o ALG8_pretransfer_H40d_DSG_eq1.out -r ALG8_pretransfer_H40d_DSG_eq1.rst -ref ./ALG8_pretransfer_H40d_DSG_min2.rst -x eq1.nc
pmemd.cuda -O -i ./input_files/eq2.in -p ./ALG8_pretransfer_H40d_DSG.prmtop -c ALG8_pretransfer_H40d_DSG_eq1.rst -o ALG8_pretransfer_H40d_DSG_eq2.out -r ALG8_pretransfer_H40d_DSG_eq2.rst -ref ALG8_pretransfer_H40d_DSG_eq1.rst -x eq2.nc
pmemd.cuda -O -i ./input_files/eq3.in -p ./ALG8_pretransfer_H40d_DSG.prmtop -c ALG8_pretransfer_H40d_DSG_eq2.rst -o ALG8_pretransfer_H40d_DSG_eq3.out -r ALG8_pretransfer_H40d_DSG_eq3.rst -ref ALG8_pretransfer_H40d_DSG_eq2.rst -x eq3.nc
pmemd.cuda -O -i ./input_files/eq4.in -p ./ALG8_pretransfer_H40d_DSG.prmtop -c ALG8_pretransfer_H40d_DSG_eq3.rst -o ALG8_pretransfer_H40d_DSG_eq4.out -r ALG8_pretransfer_H40d_DSG_eq4.rst -ref ALG8_pretransfer_H40d_DSG_eq3.rst -x eq4.nc
pmemd.cuda -O -i ./input_files/eq5.in -p ./ALG8_pretransfer_H40d_DSG.prmtop -c ALG8_pretransfer_H40d_DSG_eq4.rst -o ALG8_pretransfer_H40d_DSG_eq5.out -r ALG8_pretransfer_H40d_DSG_eq5.rst -ref ALG8_pretransfer_H40d_DSG_eq4.rst -x eq5.nc
pmemd.cuda -O -i ./input_files/eq6.in -p ./ALG8_pretransfer_H40d_DSG.prmtop -c ALG8_pretransfer_H40d_DSG_eq5.rst -o ALG8_pretransfer_H40d_DSG_eq6.out -r ALG8_pretransfer_H40d_DSG_eq6.rst -ref ALG8_pretransfer_H40d_DSG_eq5.rst -x eq6.nc
# 1st production
pmemd.cuda -O -i ./input_files/mdh.in -p ./ALG8_pretransfer_H40d_DSG.prmtop -c ALG8_pretransfer_H40d_DSG_eq6.rst -o ALG8_pretransfer_H40d_DSG_prod1.out -r ALG8_pretransfer_H40d_DSG_prod1.rst -ref ALG8_pretransfer_H40d_DSG_eq6.rst -x ALG8_pretransfer_H40d_DSG_prod1.nc


