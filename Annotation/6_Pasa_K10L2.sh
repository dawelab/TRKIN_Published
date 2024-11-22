#!/bin/bash
#SBATCH --job-name=Pasa_K10L2
#SBATCH --output=Pasa_K10L2_try5.out
#SBATCH --partition=batch_30d
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200gb
#SBATCH --time=240:00:00

#Load the modules 
module load PASA/2.5.3-foss-2022a
module load AGAT/1.1.0
module load Perl/5.34.1-GCCcore-11.3.0

#Define the variables
OUTDIR=/scratch/mjb51923/TRKIN_CRISPR/out_paper
#mkdir $OUTDIR/Pasa_K10L2
cd $OUTDIR/Pasa_K10L2
REF=$OUTDIR/RepeatMasker/CI66_K10L2_v1.fasta.masked

#This loads the module again, I think there are conflicting dependencies so these need to be run immidiatly before the command 
#module load AGAT/1.1.0
#Convert the Stringtie gtf to a gff file 
#agat_convert_sp_gxf2gxf.pl -g $OUTDIR/Stringtie_K10L2/K10L2_Stringtie_Assembly.gtf -o $OUTDIR/Stringtie_K10L2/K10L2_Stringtie_Assembly.gff3
#Extract the mRNA sequences from the stringtie annotation file 
#agat_sp_extract_sequences.pl -g $OUTDIR/Stringtie_K10L2/K10L2_Stringtie_Assembly.gff3 -f $REF -t exon --merge -o $OUTDIR/Stringtie_K10L2/K10L2_Stringtie_Assembly.fasta

#Combine the two transcriptome assemblies
#cat $OUTDIR/Trinity_Denovo_K10L2/trinity.K10L2.AllReps.Trinity.fasta $OUTDIR/Stringtie_K10L2/K10L2_Stringtie_Assembly.fasta > $OUTDIR/Pasa_K10L2/transcripts_K10L2.fasta

#This loads the module again, I think there are conflicting dependencies so these need to be run immidiatly before the command 
#module load PASA/2.5.3-foss-2022a
#module load Perl/5.34.1-GCCcore-11.3.0

#This creates a list of the names of the denovo assembled transcripts 
#This is the path to the installed pasa module you can find this by loading the module and running which pasa
PASA_HOME=/apps/eb/PASA/2.5.3-foss-2022a
#$PASA_HOME/misc_utilities/accession_extractor.pl < $OUTDIR/Trinity_Denovo_K10L2/trinity.K10L2.AllReps.Trinity.fasta > $OUTDIR/Pasa_K10L2/tdn.accs

#This loads the module again, I think there are conflicting dependencies so these need to be run immidiatly before the command 
#module load AGAT/1.1.0
#Convert the K10L2 braker output to a gff3 file
#agat_convert_sp_gxf2gxf.pl -g $OUTDIR/K10L2_Braker/braker.gtf -o $OUTDIR/K10L2_Braker/braker.gff3

#This loads the module again, I think there are conflicting dependencies so these need to be run immidiatly before the command 
module load PASA/2.5.3-foss-2022a
module load Perl/5.34.1-GCCcore-11.3.0
#Run the initial pasa step which aligns the transcripts to the reference., -C creates the db, -R runs the pipeline, -t is the transcript file, --TDN is the full length transcripts from the denovo assembly, --aligners is the aligners to use, --CPU is the threads this cannot be changed or it throws errors. 
#Launch_PASA_pipeline.pl -c $OUTDIR/Pasa_K10L2/alignAssembly.config -C -R -g $REF -t $OUTDIR/Pasa_K10L2/transcripts_K10L2.fasta --TDN $OUTDIR/Pasa_K10L2/tdn.accs --ALIGNERS minimap2 --CPU 2

#Generate the comprehensive transcriptome database
#$PASA_HOME/scripts/build_comprehensive_transcriptome.dbi -c $OUTDIR/Pasa_K10L2/alignAssembly.config -t $OUTDIR/Pasa_K10L2/transcripts_K10L2.fasta --min_per_ID 95 --min_per_aligned 30

#Load the braker genome annotations
#$PASA_HOME/scripts/Load_Current_Gene_Annotations.dbi -c $OUTDIR/Pasa_K10L2/alignAssembly.config -g $REF -P $OUTDIR/K10L2_Braker/braker.gff3

#Validate the gff3 file 
#$PASA_HOME/misc_utilities/pasa_gff3_validator.pl $OUTDIR/K10L2_Braker/braker.gff3

#Update the annotations with the transcript data
#Launch_PASA_pipeline.pl -c $OUTDIR/Pasa_K10L2/annotationCompare.config -A -g $REF -t $OUTDIR/Pasa_K10L2/transcripts_K10L2.fasta --CPU 24

#Load the updated genome anotation
#$PASA_HOME/scripts/Load_Current_Gene_Annotations.dbi -c $OUTDIR/Pasa_K10L2/alignAssembly.config -g $REF -P $OUTDIR/Pasa_K10L2/pasa_db_K10L2.gene_structures_post_PASA_updates.3434867.gff3

#Update the annotations with the transcript data again
#Launch_PASA_pipeline.pl -c $OUTDIR/Pasa_K10L2/annotationCompare.config -A -g $REF -t $OUTDIR/Pasa_K10L2/transcripts_K10L2.fasta --CPU 24

#Load the updated genome anotation
$PASA_HOME/scripts/Load_Current_Gene_Annotations.dbi -c $OUTDIR/Pasa_K10L2/alignAssembly.config -g $REF -P $OUTDIR/Pasa_K10L2/CI66_K10L2.genes.gff3

#Update the annotations with the transcript data again
Launch_PASA_pipeline.pl -c $OUTDIR/Pasa_K10L2/annotationCompare.config -A -g $REF -t $OUTDIR/Pasa_K10L2/transcripts_K10L2.fasta --CPU 24