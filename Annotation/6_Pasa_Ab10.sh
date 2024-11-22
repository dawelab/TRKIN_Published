#!/bin/bash
#SBATCH --job-name=Pasa_Ab10
#SBATCH --output=Pasa_Ab10_try4.out
#SBATCH --partition=highmem_30d_p
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=400gb
#SBATCH --time=288:00:00

#Load the modules 
module load PASA/2.5.3-foss-2022a
module load AGAT/1.1.0
module load Perl/5.34.1-GCCcore-11.3.0

#Define the variables
OUTDIR=/scratch/mjb51923/TRKIN_CRISPR/out_paper
#mkdir $OUTDIR/Pasa_Ab10
cd $OUTDIR/Pasa_Ab10
REF=$OUTDIR/RepeatMasker/Ab10_HiFi_v2_corrected.fa.masked

#This loads the module again, I think there are conflicting dependencies so these need to be run immidiatly before the command 
#module load AGAT/1.1.0
#Convert the Stringtie gtf to a gff file 
#agat_convert_sp_gxf2gxf.pl -g $OUTDIR/Stringtie_Ab10/B73Ab10_Stringtie_Assembly.gtf -o $OUTDIR/Stringtie_Ab10/B73Ab10_Stringtie_Assembly.gff3
#Extract the mRNA sequences from the stringtie annotation file 
#agat_sp_extract_sequences.pl -g $OUTDIR/Stringtie_Ab10/B73Ab10_Stringtie_Assembly.gff3 -f $REF -t exon --merge -o $OUTDIR/Stringtie_Ab10/B73Ab10_Stringtie_Assembly.fasta

#Combine the two transcriptome assemblies
#cat $OUTDIR/Trinity_Denovo_Ab10/trinity.B73Ab10.AllTissues.Trinity.fasta $OUTDIR/Stringtie_Ab10/B73Ab10_Stringtie_Assembly.fasta > $OUTDIR/Pasa_Ab10/transcripts_B73Ab10.fasta

#This loads the module again, I think there are conflicting dependencies so these need to be run immidiatly before the command 
#module load PASA/2.5.3-foss-2022a
#module load Perl/5.34.1-GCCcore-11.3.0
#This creates a list of the names of the denovo assembled transcripts 
#This is the path to the installed pasa module you can find this by loading the module and running which pasa
PASA_HOME=/apps/eb/PASA/2.5.3-foss-2022a
#$PASA_HOME/misc_utilities/accession_extractor.pl < $OUTDIR/Trinity_Denovo_Ab10/trinity.B73Ab10.AllTissues.Trinity.fasta > $OUTDIR/Pasa_Ab10/tdn.accs

#This loads the module again, I think there are conflicting dependencies so these need to be run immidiatly before the command 
#module load AGAT/1.1.0
#Convert the Ab10 braker output to a gff3 file
#agat_convert_sp_gxf2gxf.pl -g $OUTDIR/Braker/braker.gtf -o $OUTDIR/Braker/braker.gff3

#This loads the module again, I think there are conflicting dependencies so these need to be run immidiatly before the command 
module load PASA/2.5.3-foss-2022a
module load Perl/5.34.1-GCCcore-11.3.0
#Run the initial pasa step which aligns the transcripts to the reference., -C creates the db, -R runs the pipeline, -t is the transcript file, --TDN is the full length transcripts from the denovo assembly, --aligners is the aligners to use, --CPU is the threads this cannot be changed or it throws errors. 
#Launch_PASA_pipeline.pl -c $OUTDIR/Pasa_Ab10/alignAssembly.config -C -R -g $REF -t $OUTDIR/Pasa_Ab10/transcripts_B73Ab10.fasta --TDN $OUTDIR/Pasa_Ab10/tdn.accs --ALIGNERS minimap2 --CPU 2

#Generate the comprehensive transcriptome database
#$PASA_HOME/scripts/build_comprehensive_transcriptome.dbi -c $OUTDIR/Pasa_Ab10/alignAssembly.config -t $OUTDIR/Pasa_Ab10/transcripts_B73Ab10.fasta --min_per_ID 95 --min_per_aligned 30

#Load the braker genome annotations
#$PASA_HOME/scripts/Load_Current_Gene_Annotations.dbi -c $OUTDIR/Pasa_Ab10/alignAssembly.config -g $REF -P $OUTDIR/Braker/braker.gff3

#Validate the gff3 file 
#$PASA_HOME/misc_utilities/pasa_gff3_validator.pl $OUTDIR/Braker/braker.gff3

#Update the annotations with the transcript data
#Launch_PASA_pipeline.pl -c $OUTDIR/Pasa_Ab10/annotationCompare.config -A -g $REF -t $OUTDIR/Pasa_Ab10/transcripts_B73Ab10.fasta --CPU 24

#Load the updated genome anotation
#$PASA_HOME/scripts/Load_Current_Gene_Annotations.dbi -c $OUTDIR/Pasa_Ab10/alignAssembly.config -g $REF -P $OUTDIR/Pasa_Ab10/pasa_db_Ab10.gene_structures_post_PASA_updates.2685167.gff3

#Update the annotations with the transcript data again
#Launch_PASA_pipeline.pl -c $OUTDIR/Pasa_Ab10/annotationCompare.config -A -g $REF -t $OUTDIR/Pasa_Ab10/transcripts_B73Ab10.fasta --CPU 24

#Load the updated genome anotation
$PASA_HOME/scripts/Load_Current_Gene_Annotations.dbi -c $OUTDIR/Pasa_Ab10/alignAssembly.config -g $REF -P $OUTDIR/Pasa_Ab10/HiFiAb10.genes.gff3

#Update the annotations with the transcript data again
Launch_PASA_pipeline.pl -c $OUTDIR/Pasa_Ab10/annotationCompare.config -A -g $REF -t $OUTDIR/Pasa_Ab10/transcripts_B73Ab10.fasta --CPU 24


