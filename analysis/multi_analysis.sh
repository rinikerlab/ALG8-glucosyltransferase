# this script performs multiple analyses on a given dataset
#!/bin/bash


MD_path_pretransfer="/localhome/shuchen/projects/ALG8/pre_transfer/stripped/"
MD_path_posttransfer="/localhome/shuchen/projects/ALG8/post_transfer/simulations/stripped/"
H40_state=("none" "H40d" "H40db")
N_run=5
### Pre-transfer analyses
for H40 in "${H40_state[@]}"; do
    for sub in "DSG" "DSM"; do
        if [ "$H40" == "none" ]; then
            prefix="ALG8_pretransfer_${sub}"
        else
            prefix="ALG8_pretransfer_${H40}_${sub}"
        fi
        python ./ALG8_analysis.py --MD_path $MD_path_pretransfer --prefix $prefix --N_run $N_run --state "pre-$sub" --take_frame 4000
        # echo "python ./ALG8_analysis.py --MD_path $MD_path_pretransfer --prefix $prefix --N_run $N_run --state pre-$sub --take_frame 4000"
    done
done

### Post-transfer analyses
for H40 in "${H40_state[@]}"; do
    for act in "D36" "N36"; do
        if [ "$H40" == "none" ]; then
            prefix="ALG8_${act}_post_transfer"
        else
            prefix="ALG8_${act}_${H40}_post_transfer"
        fi
        python ./ALG8_analysis.py --MD_path $MD_path_posttransfer --prefix $prefix --N_run $N_run --state "post" --take_frame 4000
    done
done

# DSG
