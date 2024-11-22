#!/bin/bash
#SBATCH --job-name=EnTAP_Ab10
#SBATCH --output=EnTAP_Ab10.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200gb
#SBATCH --time=100:00:00

#load the modules
module load EnTAP/1.0.0-foss-2022a
module load SAMtools/1.14-GCC-11.2.0

#Define variables
DIR=/scratch/mjb51923/TRKIN_CRISPR/out_paper/Entap

#Do the basic EnTAP configuration
cd $DIR/Ab10/
#EnTAP --config --ini /home/mjb51923/entap_config.ini

#Configure my databases to be diamond databases
######################you will need to redo this when the nr database is ready 
#To get EnTAP running you need to download these databases and then run the following commands. Refer to here for instructions https://entap.readthedocs.io/en/latest/Getting_Started/Configuration/configuration.html

EnTAP --config --ini /home/mjb51923/entap_config.ini -i /scratch/mjb51923/TRKIN_CRISPR/out_paper/AGAT/HiFiAb10.v2.cds.fasta -d /scratch/mjb51923/TRKIN_CRISPR/out_paper/Entap/uniprot_sprot.fasta -d /scratch/mjb51923/TRKIN_CRISPR/out_paper/Entap/RefSeq.plant.protein.faa -d /scratch/mjb51923/TRKIN_CRISPR/out_paper/Entap/nr.fasta --overwrite

#This runs entap using the diamond databases from the above step. runP means it will do frame selection, convert to protein, and then use blastp against the databases as I am providing a CDS transcriptome
module load EnTAP/1.0.0-foss-2022a
EnTAP --runP --ini /home/mjb51923/entap_config.ini -i /scratch/mjb51923/TRKIN_CRISPR/out_paper/AGAT/HiFiAb10.v2.cds.fasta -d /scratch/mjb51923/TRKIN_CRISPR/out_paper/Entap/K10L2/entap_outfiles/bin/uniprot_sprot.dmnd -d /scratch/mjb51923/TRKIN_CRISPR/out_paper/Entap/K10L2/entap_outfiles/bin/RefSeq.dmnd -d /scratch/mjb51923/TRKIN_CRISPR/out_paper/Entap/K10L2/entap_outfiles/bin/nr.dmnd  --out-dir /scratch/mjb51923/TRKIN_CRISPR/out_paper/Entap/Ab10/entap_outfiles_v2 -t 24
