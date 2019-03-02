#!/bin/bash

pwd

#output to PED
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --recode --tab --out ../MIPS_SexCorrected

#/storage/software/plink --file /storage/mips/MIPS_Updated.2019-02-21/data/ --pheno 20190301_tcrEmrPheno.txt

/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --pheno ../20190301_tcrEmrPheno.txt --pheno-name Productive.Clonality -assoc --make-bed --out ../MIPS_SexCorrected_Pheno

