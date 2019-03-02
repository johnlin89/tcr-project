#!/bin/bash

pwd
outputFolder="/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/"

# Output PED
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --recode --tab --out $outputFolder/MIPS_SexCorrectedPed

# Add in phenotypes
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --pheno ../20190301_tcrEmrPheno.txt --pheno-name Productive.Clonality -assoc --make-bed --out $outputFolder/phenotypeTest/MIPS_SexCorrected_Pheno

# Linear Model with TCRB filtering
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --pheno ../20190301_tcrEmrPheno.txt --pheno-name Productive.Clonality --chr 7 --from-mb 141.2 --to-kb 142 --linear --make-bed --out $outputFolder/linearTest/MIPS_SexCorrected_Pheno_tcrb_linear
