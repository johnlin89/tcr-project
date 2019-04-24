#!/bin/bash

pwd
outputFolder="/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/"

# updated coordinates for TRB locus, also updated for maf and hwe filtering all in one shot 
rm -r $outputFolder/plinkFiltering/plink1
mkdir $outputFolder/plinkFiltering/plink1
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --pheno ../phenotypes/tcrEmrPheno.txt --pheno-name Productive.Clonality --chr 7 --from-bp 141998851 --to-bp 142510972 --maf 0.1 --hwe 0.0001 --linear --make-bed --out $outputFolder/plinkFiltering/plink1/updatedCoord
# 600 SNPs in region (1824820 snps total)
# 489 removed due to MAF, 3 due to HWE

# just output after MAF and HWE Filtering
rm -r $outputFolder/plinkFiltering/plink2
mkdir $outputFolder/plinkFiltering/plink2
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --pheno ../phenotypes/tcrEmrPheno.txt --pheno-name Productive.Clonality --chr 7 --from-bp 141998851 --to-bp 142510972 --maf 0.1 --hwe 0.0001 --make-bed --out $outputFolder/plinkFiltering/plink2/newFilter --freq counts
# 600 SNPs in region 
# 489 removed due to MAF, 3 due to HWE 
# leaving 108 variants

# prune
rm -r $outputFolder/plinkFiltering/plink3
mkdir $outputFolder/plinkFiltering/plink3
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --pheno ../phenotypes/tcrEmrPheno.txt --pheno-name Productive.Clonality --chr 7 --from-bp 141998851 --to-bp 142510972 --maf 0.1 --hwe 0.0001 --make-bed --out $outputFolder/plinkFiltering/plink3/prune --freq counts --prune
# 600 SNPs in region 
# 459 removed due to MAF, 0 due to HWE 
# leaving 141 variants

# analysis: do linear regressions after pruning for 15
rm -r $outputFolder/plinkFiltering/plink4
mkdir $outputFolder/plinkFiltering/plink4
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink3/prune --pheno ../phenotypes/tcrEmrPheno.txt --pheno-name Productive.Clonality --chr 7 --from-bp 141998851 --to-bp 142510972 --linear --make-bed --out $outputFolder/plinkFiltering/plink4/pruneLinear --freq counts
# 141 variants

# expand window by 50000 bp on each side and do linear regressions
rm -r $outputFolder/plinkFiltering/plink5
mkdir $outputFolder/plinkFiltering/plink5
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --pheno ../phenotypes/tcrEmrPheno.txt --pheno-name Productive.Clonality --chr 7 --from-bp 141948851 --to-bp 142560972  --maf 0.1 --hwe 0.0001 --linear --make-bed --out $outputFolder/plinkFiltering/plink5/window --freq counts --prune
# 654 snps in window
# 495 removed due to MAF, 1 removed due to HWE
# 159 SNPs left


# see all of chr 7 with pruned individuals, and MAF and HWE Filtering
rm -r $outputFolder/plinkFiltering/plink6
mkdir $outputFolder/plinkFiltering/plink6
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --pheno ../phenotypes/tcrEmrPheno.txt --pheno-name Productive.Clonality --chr 7 --maf 0.1 --hwe 0.0001 --linear --make-bed --out $outputFolder/plinkFiltering/plink6/allChr7 --freq counts --prune
# 97332 snps in window
# 63571 removed due to MAF, 1 removed due to HWE
# 33760 SNPs left

# use age and sex as covariates and redo linear regressions, building off of expanded window from plink5
rm -r $outputFolder/plinkFiltering/plink7
mkdir $outputFolder/plinkFiltering/plink7
/storage/software/plink --bfile $outputFolder/plinkFiltering/plink5/window --covar ../phenotypes/tcrEmrPheno.txt  --covar-name SEX,AGE --linear --make-bed --out $outputFolder/plinkFiltering/plink7/windowCovar
# 159 SNPs

# use all samples (not just 15) to estimate admixture
# use expanded window, 50,000 bp on each side of TRB locus
# got rid of pheno options and linear modeling
rm -r $outputFolder/plinkFiltering/plink8
mkdir $outputFolder/plinkFiltering/plink8
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --chr 7 --from-bp 141948851 --to-bp 142560972  --maf 0.1 --hwe 0.0001 --make-bed --out $outputFolder/plinkFiltering/plink8/allWindow --freq counts
# 129 people (61 males, 68 females)
# 654 snps from file
# 3 snps removed due to hwe
# 527 removed due to MAF
# 124 variants removed