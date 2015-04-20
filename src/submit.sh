#!/bin/bash

#$ -N find_phase
#$ -o find_phase.out
#$ -e find_phase.error
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=156:00:00,h_vmem=5G


python do_it.py