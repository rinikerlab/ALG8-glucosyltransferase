#!/bin/bash
# Copyright (C) 2026 ETH Zurich, Shu-Yu Chen
# this script performs multiple analyses on a given dataset
stripped_home='../data/'

N_run=5

# pretransfer
folders=( ALG8_pretransfer_H40d_DSG ALG8_pretransfer_H40e_DSG ALG8_pretransfer_H40d_DSM ALG8_pretransfer_H40e_DSM)
prefix=( ALG8_pretransfer_H40d_DSG ALG8_pretransfer_DSG ALG8_pretransfer_H40d_DSM ALG8_pretransfer_DSM)
for i in $(seq 0 3); do
    # python ./ALG8_analysis.py --MD_path "${stripped_home}/${folders[i]}/" --prefix "${prefix[i]}" --N_run $N_run --state "pre" --take_frame 6000
    echo "python ./ALG8_analysis.py --MD_path ${stripped_home}/${folders[i]}/ --prefix ${prefix[i]} --N_run $N_run --state pre --take_frame 4000"
done

# postransfer
folders=( ALG8_postransfer_H40d ALG8_postransfer_H40e  ALG8_postransfer_H40d_D36N ALG8_postransfer_H40e_D36N)
prefix=( ALG8_D36_H40d_post_transfer ALG8_D36_post_transfer ALG8_N36_H40d_post_transfer ALG8_N36_post_transfer)
for i in $(seq 0 3); do
    python ./ALG8_analysis.py --MD_path "${stripped_home}/${folders[i]}/" --prefix "${prefix[i]}" --N_run $N_run --state "post" --take_frame 7500
    # echo "python ./ALG8_analysis.py --MD_path ${stripped_home}/${folders[i]}/ --prefix ${prefix[i]} --N_run $N_run --state post --take_frame 4000"
done

