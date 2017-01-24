#!/bin/bash 

#SBATCH --time=00:10:00 
#SBATCH --partition=compute 
#SBATCH --mem-per-cpu=1G 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=2 
#SBATCH --job-name=sBCBG_2017_1_24_13:12_00000
#SBATCH --input=none
#SBATCH --output="2017_1_24_13:12_00000.out" 
#SBATCH --error="2017_1_24_13:12_00000.err" 
#SBATCH --mail-user=benoit.girard@isir.upmc.fr 
#SBATCH --mail-type=BEGIN,END,FAIL 
module load nest/2.10 
python testFullBG.py 
