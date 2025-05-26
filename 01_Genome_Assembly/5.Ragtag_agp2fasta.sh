#load the modules
module load RagTag/2.1.0

OUTDIR=""

#Assemble the K10L2 while leaving all other contigs the same. I determined this order by scaffolding on B73 and verifying the position of K10L2 landmarks.
ragtag.py agp2fa $OUTDIR/ragtag.scaffold.K10L2.agp $OUTDIR/CI66_rq99.asm.bp.p_ctg.gfa.fasta > $OUTDIR/RepeatMasker/CI66_K10L2_v1.fasta
