module load AGAT/1.1.0
#Convert the gtf file from Liftoff to a gff file. This file is from step 05.08
gunzip /path/to/Liftoff/B73_Ab10_HiFi_v2.gene.gff3.gz
agat_convert_sp_gff2gtf.pl --gff /path/to/Liftoff/B73_Ab10_HiFi_v2.gene.gff3 -o /path/to/Liftoff/B73_Ab10_HiFi_v2.gene.gtf
