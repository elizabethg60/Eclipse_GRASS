#!/bin/bash
#SBATCH --account=sqm107_b
#SBATCH --partition=sla-prio
#SBATCH --time=48:00:00
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --mem=10gb
#SBATCH --job-name=figures
#SBATCH --output /storage/home/efg5335/work/Eclipse_GRASS/figures.out
#SBATCH --error /storage/home/efg5335/work/Eclipse_GRASS/figures.err
#SBATCH --mail-user=efg5335@psu.com

##now=$(date +"%T")
##echo "start time ... $now"

cd /storage/home/efg5335/work/Eclipse_GRASS
julia time_loop.jl

##now=$(date +"%T")
##echo "end time ... $now"