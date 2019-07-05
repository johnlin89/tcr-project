#!/bin/bash

# For Pacific Symposium on Biocomputing Precision Medicine Category
# This is to identify the SNPs in the MIPs dataset that are found to have some clinical significance from ClinVar

pwd
outputFolder = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults"
referenceFolder = "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/PAGEII_WGSA_MEGA_Annotations/mega_snp"

# Use plink to output all snps with just genotyping rate and HWE filters
rm -r $outputFolder/plinkFiltering/plink16
mkdir $outputFolder/plinkFiltering/plink16
/storage/software/plink --bfile /storage/mips/MIPS_Updated.2019-02-21/data/MIPS_SexCorrected --geno 0.02 --hwe 0.0001 --make-bed --out $outputFolder/plinkFiltering/plink16/clinVar
# 1790770 variants left

# how many snps in the sorted_mega_annotated.snp file?
echo sorted_mega_annotated.snp has this many lines:
wc -l sorted_mega_annotated.snp
# 1372956 records in sorted_mega_annotated.snp

# From the clinvar data - filter for just those records that have clinvar_rs column populated
echo Filtering sorted_mega_annotated.snp for populated clinvar_rs column...
awk -F "\t" '{ if ($139 != ".") { print } }' $referenceFolder/sorted_mega_annotated.snp | cut -f 1,2,3,4,5,6,139,140,141,142 > $referenceFolder/sorted_mega_annotated_clinVar.snp
echo After filtering, sorted_mega_annotated_clinVar.snp has this many lines:
wc -l sorted_mega_annotated_clinVar.snp
# 23720 snps left after filtering for clinvar_rs

# sort by snp id so that we can join later
echo Formatting and sorting clinVar.bim for join...
echo -e "chr\trsid\tgd\tbp\ta1\ta2" | cat - $outputFolder/plinkFiltering/plink16/clinVar.bim > $outputFolder/plinkFiltering/plink16/clinVarHeader.bim
head -n 1 $outputFolder/plinkFiltering/plink16/clinVarHeader.bim > $outputFolder/plinkFiltering/plink16/clinVarSort.bim &&
tail -n +2 $outputFolder/plinkFiltering/plink16/clinVarHeader.bim | sort -k 2 >> $outputFolder/plinkFiltering/plink16/clinVarSort.bim
echo Sorting sorted_mega_annotated_clinVar for join...
head -n 1 $referenceFolder/sorted_mega_annotated_clinVar.snp > $referenceFolder/sorted_mega_annotated_clinVar_sorted.snp &&
tail -n +2 $referenceFolder/sorted_mega_annotated_clinVar.snp | sort -k 7 >> $referenceFolder/sorted_mega_annotated_clinVar_sorted.snp

echo Joining clinVarSort.bim and sorted_mega_annotated_clinVar_sorted.snp
join -1 2 -2 5 $outputFolder/plinkFiltering/plink16/clinVarSort.bim $referenceFolder/sorted_mega_annotated_clinVar_sorted.snp > $outputFolder/mipsClinVar.txt
echo This many snps were found between the PAGEII_WGSA_MEGA_Annotations and the MIPS genotype data:
wc -l mipsClinVar.txt
# 1179 snps
# which people have which genotypes??