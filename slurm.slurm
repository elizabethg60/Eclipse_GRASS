#!/bin/bash
#SBATCH --account=sqm107_b
#SBATCH --time=72:00:00
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --mem=15gb
#SBATCH --job-name=my_model_full
#SBATCH --output /storage/home/efg5335/work/Eclipse_GRASS/my_model_full.out
#SBATCH --error /storage/home/efg5335/work/Eclipse_GRASS/my_model_full.err
#SBATCH --mail-user=efg5335@psu.com

##now=$(date +"%T")
##echo "start time ... $now"

cd /storage/home/efg5335/work/Eclipse_GRASS
julia my_model.jl

##now=$(date +"%T")
##echo "end time ... $now"