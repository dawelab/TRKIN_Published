#Load the modules
module load Trinity/2.15.1-foss-2022a

#Define the variables
OUTDIR=""
mkdir $OUTDIR/Trinity_Denovo_K10L2
READDIR=/path/to/RNAseq_K10L2

#List the fastq files 
R1_LIST=$(ls $READDIR/*_R1_paired.fq)
R2_LIST=$(ls $READDIR/*_R2_paired.fq)

#Concatenate the fastq reads
cat $R1_LIST > $OUTDIR/Trinity_Denovo_K10L2/K10L2.AllReps_R1.fq
cat $R2_LIST > $OUTDIR/Trinity_Denovo_K10L2/K10L2.AllReps_R2.fq

#Denovo assemble the K10L2 reads
for i in "K10L2.AllReps"
do
echo "########### Starting new line"
mkdir $OUTDIR/Trinity_Denovo_K10L2/trinity.$i
#Perform denovo transcriptome assembly using all available RNAseq data
Trinity --seqType fq --output $OUTDIR/Trinity_Denovo_K10L2/trinity.$i --max_memory 170G --left $OUTDIR/Trinity_Denovo_K10L2/${i}_R1.fq --right $OUTDIR/Trinity_Denovo_K10L2/${i}_R2.fq --CPU 24
done
