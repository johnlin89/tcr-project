referenceFolder="/storage/mips/MIPS_Updated.2019-02-21/data/1KG_Phase3"
outputFolder = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults"

# From CEU-YRI_subsetter.sh
# Preparation of 1000 Genomes for merger with local files 

### You may need to hardcode plink2's location if not in your PATH and modify reference directory ###
#module load plink2
ref_dir="/storage/mips/MIPS_Updated.2019-02-21/data/1KG_Phase3"

# Create list of populations and subset 
# awk '{if($6=="CEU" || $6=="YRI") print}' $ref_dir/all_phase3.psam > $ref_dir/CEU-YRI.psam 

# add #IID	PAT	MAT	SEX	SuperPop	Population
# to first line of CEU-YRI.psam manually

# Omit multiallelic, limit to populations of interest, create bed 
# /storage/software/plink2 -pfile $ref_dir/all_phase3 \
#  --chr 1-22 \
#  --keep $ref_dir/CEU-YRI.psam \
#  --make-bfile --out $ref_dir/temp

# 2504 samples (1271 females, 1233 males; 2497 founders) loaded from
# /storage/mips/MIPS_Updated.2019-02-21/data/1KG_Phase3/all_phase3.psam.
# 81271745 out of 84805772 variants loaded from
# /storage/mips/MIPS_Updated.2019-02-21/data/1KG_Phase3/all_phase3.pvar.
# 207 samples (106 females, 101 males; 207 founders) remaining after main filters.

# Removal of 1st and 2nd degree related samples using KING file from Plink2's website
# /storage/software/plink2 -bfile $ref_dir/temp --remove $ref_dir/deg2.king.cutoff.out.id --make-bed --out $ref_dir/CEU-YRI_unrel

# no people are removed after king relatedness above

# Cleanup
# rm $ref_dir/temp*

# filter for chr 1-22
# bi allelic only
# /storage/software/plink --bfile $ref_dir/CEU-YRI_unrel --chr 1-22 --biallelic-only strict --make-bed --out $referenceFolder/CEU-YRI_unrelFinal

# exclude some snps
/storage/software/plink -bfile $referenceFolder/CEU-YRI_unrelFinal --exclude /storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink10/allWindowNoDup_flipRef-merge.missnp --make-bed --out $referenceFolder/CEU-YRI_unrelFinal2