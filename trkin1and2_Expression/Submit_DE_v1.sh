#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Submit_DE_v1
#SBATCH --output Submit_DE_v1.out
#SBATCH --mem=100GB
#SBATCH --time=24:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=12
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#Load the modules needed for sapelo2
module load R/4.3.1-foss-2022a

#set the working directory 
cd /scratch/mjb51923/TRKIN_CRISPR/TRKIN_Published/trkin1and2_Expression

R CMD BATCH DE_v1.R
