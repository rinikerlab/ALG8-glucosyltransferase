# The scripts here are shortened by 100 times for demo. For original simulation length please replace the input files with the ones in the longer_intput_files folder.

# minimization
bash min.sh

# equilibration
bash eq.sh

# production for 5ns (1ns * 5) 
for i in $(seq 1 4)
do
bash pro.sh $i
done
