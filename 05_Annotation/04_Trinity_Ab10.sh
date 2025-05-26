#Load the modules
module load Trinity/2.15.1-foss-2022a

#Define the variables
OUTDIR=""
#mkdir $OUTDIR/Trinity_Denovo_Ab10
READDIR=/path/to/RNAseq_Ab10

#List the fastq files 
R1_LIST=$(ls $READDIR/*_R1_paired.fq)
R2_LIST=$(ls $READDIR/*_R2_paired.fq)

#Concatenate the fastq reads
cat $R1_LIST > $OUTDIR/Trinity_Denovo_Ab10/B73Ab10.AllTissues_R1.fq
cat $R2_LIST > $OUTDIR/Trinity_Denovo_Ab10/B73Ab10.AllTissues_R2.fq

#Denovo assemble the Ab10 reads
for i in "B73Ab10.AllTissues"
do
echo "########### Starting new line"
mkdir $OUTDIR/Trinity_Denovo_Ab10/trinity.$i
#This performs the denovo transcriptome assembly using all available RNAseq data
Trinity --seqType fq --output $OUTDIR/Trinity_Denovo_Ab10/trinity.$i --max_memory 270G --left $OUTDIR/Trinity_Denovo_Ab10/${i}_R1.fq --right $OUTDIR/Trinity_Denovo_Ab10/${i}_R2.fq --CPU 24
done
