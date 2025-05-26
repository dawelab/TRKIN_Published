#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Ragtag_agp2fasta
#SBATCH --output Ragtag_agp2fasta.out
#SBATCH --mem=2GB
#SBATCH --time=5:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#load the modules
module load RagTag/2.1.0

#Assemble the K10L2 while leaving all other contigs the same. I determined this order by scaffolding on B73 and verifying the position of K10L2 landmarks.
ragtag.py agp2fa /scratch/mjb51923/TRKIN_CRISPR/TRKIN_Paper/ragtag.scaffold.K10L2.agp /scratch/mjb51923/CI66_Assembly/out/CI66_rq99.asm.bp.p_ctg.gfa.fasta > /scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/CI66_K10L2_v1.fasta
