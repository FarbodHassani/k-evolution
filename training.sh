#!/bin/bash

while read c1 c2 c3 c4 c5 c6; do Omega_cdm_arr[i]=$c1 n_s_arr[i]=$c2 h_arr[i]=$c3 w_kess_arr[i]=$c4 cs2_arr[i]=$c5 A_s_arr[i++]=$c6; done < LHS_125.txt


for ((i = 0 ; i < 2 ; i++)); do

Omega_cdm=${Omega_cdm_arr[$i]}
n_s=${n_s_arr[$i]}
h=${h_arr[$i]}
w_kessence=${w_kess_arr[$i]}
cs2_kessence=$(awk -v var="${cs2_arr[$i]}" 'BEGIN {print 10^var}')
A_s=${A_s_arr[$i]}

output_dir='output_'$i''
#echo $output_dir;

#########
[ ! -d "./LHS_outputs/$output_dir" ] && mkdir -p "./LHS_outputs/$output_dir"


sed -e 's/Omega_cdm   = 0.27/Omega_cdm   = '$Omega_cdm'/g' -e 's/n_s = 0.9619/n_s = '$n_s'/g' -e 's/h           = 0.67556/h           = '$h'/g' -e 's/w_kessence = -0.9/w_kessence = '$w_kessence'/g' -e 's/cs2_kessence = 1.e-7/cs2_kessence = '$cs2_kessence'/g' -e 's/A_s = 2.215e-9/A_s = '$A_s'/g' -e 's#output path         = output/#output path         = 'LHS_outputs/$output_dir/'#g' settings_unity.ini> settings_$output_dir.ini

sed 's/settings_unity.ini/settings_'$output_dir'.ini/g' Run_kev.sh> Run_number$i.sh
#
bash Run_number$i.sh

mv settings_$output_dir.ini ./LHS_outputs/$output_dir
mv Run_number$i.sh ./LHS_outputs/$output_dir

done
