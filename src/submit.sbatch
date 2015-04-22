#!/bin/bash

#SBATCH --job-name 128_find_phase
#SBATCH --output find_phase.%j.%N.out
#SBATCH --error find_phase.%j.%N.err
#SBATCH --get-user-env
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH -n 16
#SBATCH --mem-per-cpu=64000
#SBATCH --time=24:00:00 

source ~/code/cs_env/bin/activate
python do_it.py &
python do_it2.py
