#!/bin/bash
#SBATCH --partition=batch 
#SBATCH -J Download_Gapless_RNA
#SBATCH --output Download_Gapless_RNA.out
#SBATCH --mem=10GB
#SBATCH --time=10:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=30
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=END

cd /scratch/mjb51923/raw_reads/RNA/Gapless_B73-Ab10I

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801510/B73Ab10_V11_middle_MN02041_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801514/B73Ab10_V18_tassel_MN02061_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801519/B73Ab10_R1_anther_MN02082_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801521/B73Ab10_16DAP_endosperm_MN02092_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801513/B73Ab10_V11_tip_MN02052_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801522/B73Ab10_16DAP_embryo_MN02101_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801505/B73Ab10_8DAS_root_MN02012_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801520/B73Ab10_16DAP_endosperm_MN02091_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801523/B73Ab10_16DAP_embryo_MN02102_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801518/B73Ab10_R1_anther_MN02081_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801512/B73Ab10_V11_tip_MN02051_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801517/B73Ab10_V18_ear_MN02072_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801504/B73Ab10_8DAS_root_MN02011_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801509/B73Ab10_V11_base_MN02032_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801516/B73Ab10_V18_ear_MN02071_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801507/B73Ab10_8DAS_shoot_MN02022_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801505/B73Ab10_8DAS_root_MN02012_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801508/B73Ab10_V11_base_MN02031_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801514/B73Ab10_V18_tassel_MN02061_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801516/B73Ab10_V18_ear_MN02071_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801507/B73Ab10_8DAS_shoot_MN02022_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801521/B73Ab10_16DAP_endosperm_MN02092_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801506/B73Ab10_8DAS_shoot_MN02021_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801504/B73Ab10_8DAS_root_MN02011_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801511/B73Ab10_V11_middle_MN02042_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801515/B73Ab10_V18_tassel_MN02062_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801515/B73Ab10_V18_tassel_MN02062_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801518/B73Ab10_R1_anther_MN02081_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801509/B73Ab10_V11_base_MN02032_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801511/B73Ab10_V11_middle_MN02042_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801506/B73Ab10_8DAS_shoot_MN02021_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801519/B73Ab10_R1_anther_MN02082_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801512/B73Ab10_V11_tip_MN02051_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801523/B73Ab10_16DAP_embryo_MN02102_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801520/B73Ab10_16DAP_endosperm_MN02091_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801508/B73Ab10_V11_base_MN02031_R2.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801517/B73Ab10_V18_ear_MN02072_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801510/B73Ab10_V11_middle_MN02041_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801513/B73Ab10_V11_tip_MN02052_R1.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR380/ERR3801522/B73Ab10_16DAP_embryo_MN02101_R2.fq.gz

