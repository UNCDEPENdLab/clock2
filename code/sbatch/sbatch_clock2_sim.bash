#!/bin/bash

#SBATCH -p smp
#SBATCH -N 1
#SBATCH --mem 1g
#SBATCH -n 1
#SBATCH -t 23:00:00
#SBATCH --mail-user=ayd1@pitt.edu
#SBATCH --mail-type=ALL
#SBATCH -c 1
[ -z "$sourcefilestart" ] && echo "No sourcefilestart env variable passed in" && exit 1

cd ../ #sbatch scripts live in subdirectory
R CMD BATCH --no-save --no-restore clock2_sim_crc.R clock2_${sourcefilestart}_$(date +%s).Rout