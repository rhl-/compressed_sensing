#!/bin/bash

#$ -N b2find_phase
#$ -o b2find_phase.out
#$ -e b2find_phase.error
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=156:00:00,h_vmem=5G

module load python/2.7.3.mkl
python do_it.py
