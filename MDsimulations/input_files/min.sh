mpirun -n 8 pmemd.MPI -O -i ./min.in -p ALG8_pretransfer_H40d_DSG.prmtop -c ALG8_pretransfer_H40d_DSG.inpcrd -o ALG8_pretransfer_H40d_DSG_min.out -r ALG8_pretransfer_H40d_DSG_min.rst -ref ALG8_pretransfer_H40d_DSG.inpcrd -x min1.nc
mpirun -n 8 pmemd.MPI -O -i ./min2.in -p ALG8_pretransfer_H40d_DSG.prmtop -c ALG8_pretransfer_H40d_DSG_min.rst -o ALG8_pretransfer_H40d_DSG_min2.out -r ALG8_pretransfer_H40d_DSG_min2.rst -ref ALG8_pretransfer_H40d_DSG_min.rst -x min2.n




