#!/bin/bash
rm -rf namd_res_* avg_deco_* run_namd_* slurm-*
# The igeo values
for i in {0..1000..2}; do
  cp run_namd.py "run_namd_$i.py"
  sed -i "s/igeo =.*/igeo =$i/g;s/fgeo=.*/fgeo=$((i+25))/g" run_namd_$i.py
  sed -i "s/python run.*/python run_namd_$((i)).py/g" submit.slm
  sbatch submit.slm
  echo $i
  #sleep 2
done
