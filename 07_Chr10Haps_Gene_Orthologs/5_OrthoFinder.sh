#load the modules necessary from the cluster
module load OrthoFinder/2.5.5-foss-2023a

###############################################################################################################################################################################################
############################################################                              Ab10 v K10L2                             ############################################################
###############################################################################################################################################################################################

#Define Variables 
DIR=""
ODIRNAME="Ab10_K10L2_proteomes"
PROT1=$DIR/HiFiAb10.Ab10hapProtein.LongestIsoform.fasta
PROT2=$DIR/CI66_K10L2.K10L2hapProtein.LongestIsoform.fasta

#make a directory with the proteomes of the two lines
mkdir $DIR/$ODIRNAME
cp $PROT1 $DIR/$ODIRNAME
cp $PROT2 $DIR/$ODIRNAME

#These are the files in the folder $DIR/Ab10_N10_proteomes: Zm-B73_AB10-REFERENCE-NAM-1.0_Zm00043a.1.evd.protein.Ab10HapLongest.fa Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.Longest.fa
orthofinder -t 40 -f $DIR/$ODIRNAME

###############################################################################################################################################################################################
############################################################                              Ab10 v N10                               ############################################################
###############################################################################################################################################################################################


ODIRNAME="Ab10_N10_proteomes"
PROT1=$DIR/HiFiAb10.Ab10hapProtein.LongestIsoform.fasta
PROT2=$DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.Protein.LongestIsoform.fasta

#make a directory with the proteomes of the two lines
mkdir $DIR/$ODIRNAME
cp $PROT1 $DIR/$ODIRNAME
cp $PROT2 $DIR/$ODIRNAME

#These are the files in the folder $DIR/Ab10_N10_proteomes: Zm-B73_AB10-REFERENCE-NAM-1.0_Zm00043a.1.evd.protein.Ab10HapLongest.fa Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.Longest.fa
orthofinder -t 40 -f $DIR/$ODIRNAME

###############################################################################################################################################################################################
############################################################                              K10L2 v N10                              ############################################################
###############################################################################################################################################################################################


ODIRNAME="K10L2_N10_proteomes"
PROT1=$DIR/CI66_K10L2.K10L2hapProtein.LongestIsoform.fasta
PROT2=$DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.Protein.LongestIsoform.fasta

#make a directory with the proteomes of the two lines
mkdir $DIR/$ODIRNAME
cp $PROT1 $DIR/$ODIRNAME
cp $PROT2 $DIR/$ODIRNAME

#These are the files in the folder $DIR/Ab10_N10_proteomes: Zm-B73_AB10-REFERENCE-NAM-1.0_Zm00043a.1.evd.protein.Ab10HapLongest.fa Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.Longest.fa
orthofinder -t 40 -f $DIR/$ODIRNAME
