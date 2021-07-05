#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8

cd /home/zack_ely/denovo_2

module add r/3.5.1

Rscript scde_script_denovo.R


