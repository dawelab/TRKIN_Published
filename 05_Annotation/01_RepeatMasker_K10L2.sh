#Load Modules needed 
module load RepeatMasker/4.1.5-foss-2022a

#Set the file name 
REF=CI66_K10L2_v1.fasta
OUT=""

#Make the output folder
mkdir $OUT/RepeatMasker
cd $OUT/RepeatMasker

#Download the curated maize TE library
git clone https://github.com/oushujun/MTEC.git
mv MTEC/maizeTE02052020 ./maizeTE02052020.fasta
rm -r MTEC

#This runs repeat masker as maize on the reference using blast (-e), 40 threads (-pa), the custom librarby made above (-lib), and
#outputs into the designated output directory (-dir), soft masks (-small), and outputs a TE gff file (-gff)
RepeatMasker -e ncbi -pa 40 -lib $OUT/RepeatMasker/maizeTE02052020.fasta -dir $OUT/RepeatMasker -small -xsmall -gff $REF


