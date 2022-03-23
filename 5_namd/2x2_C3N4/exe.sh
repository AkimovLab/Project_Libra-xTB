#!/bin/bash

# The igeo values
for i in {0..975..25}; do
  cp run_namd.py "run_namd_$i.py"
  sed -i "s/igeo=.*/igeo=$i/g;s/fgeo=.*/fgeo=$((i+25))/g" run_namd_$i.py
  sed -i "s/python run.*/python run_name_$((i)).py/g" submit.slm
  sbatch submit.slm
  echo $i
  sleep 10
done
