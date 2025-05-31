#!/bin/bash
#SBATCH --job-name=EnTAP_K10L2
#SBATCH --output=EnTAP_K10L2.out
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
OUTDIR=""
DIR=$OUTDIR/Entap

#Do the basic EnTAP configuration
cd $DIR/K10L2/
EnTAP --config --ini /path/to/entap_config.ini

#To get EnTAP running you need to download these databases and then run the following commands. Refer to here for instructions https://entap.readthedocs.io/en/latest/Getting_Started/Configuration/configuration.html
#Configure entap and set up diamond databases
EnTAP --config --ini /path/to/entap_config.ini -i $OUTDIR/AGAT/CI66_K10L2_v1.cds.fasta -d $OUTDIR/Entap/uniprot_sprot.fasta -d $OUTDIR/Entap/RefSeq.plant.protein.faa -d $OUTDIR/Entap/nr.fasta --overwrite

#This runs entap using the diamond databases from the above step. runP means it will do frame selection, convert to protein, and then use blastp against the databases as I am providing a CDS transcriptome
EnTAP --runP --ini /path/to/entap_config.ini -i $OUTDIR/AGAT/CI66_K10L2_v1.cds.fasta  -d $OUTDIR/Entap/K10L2/entap_outfiles/bin/uniprot_sprot.dmnd -d $OUTDIR/Entap/K10L2/entap_outfiles/bin/RefSeq.dmnd -d $OUTDIR/Entap/K10L2/entap_outfiles/bin/nr.dmnd --out-dir $OUTDIR/Entap/K10L2/entap_outfiles_nr_update -t 24
