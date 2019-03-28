#!/bin/bash

pwd
outputFolder="/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/"

# updated coordinates for TRB locus, also updated for maf and hwe filtering all in one shot 
# /storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --pheno ../tcrEmrPheno.txt --pheno-name Productive.Clonality --chr 7 --from-bp 141998851 --to-bp 142510972 --maf 0.1 --hwe 0.0001 --assoc --qt-means --linear --make-bed --out $outputFolder/plinkFiltering/plink1

# output after MAF and HWE Filtering
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --pheno ../tcrEmrPheno.txt --pheno-name Productive.Clonality --chr 7 --from-bp 141998851 --to-bp 142510972 --maf 0.1 --hwe 0.0001 --make-bed --out $outputFolder/plinkFiltering/plink2 --freq --counts

# prune
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --pheno ../tcrEmrPheno.txt --pheno-name Productive.Clonality --chr 7 --from-bp 141998851 --to-bp 142510972 --maf 0.1 --hwe 0.0001 --make-bed --out $outputFolder/plinkFiltering/plink3 --freq --counts --prune
