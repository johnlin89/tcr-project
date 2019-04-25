# for all 129 samples
# for K in 1 2 3 4 5 6 7 8 9 10
# 	do 
# 		../admixture_macosx-1.3.0/admixture --cv allWindow_flip.bed $K | tee log${K}.out
# 	done
# 	
# grep -h CV log*.out
# K = 6 shows lowest error
outputFolder="/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink15"
ref_dir="/storage/mips/MIPS_Updated.2019-02-21/data/1KG_Phase3"

# create .pop file
awk '{ if (NR!=1) {print $1 " " $6} }' $ref_dir/CEU-YRI.psam > $outputFolder/temp.pop
sort $outputFolder/temp.pop -o $outputFolder/temp.pop
sort $outputFolder/allWindowNoDup_flipRefPruned.fam -o $outputFolder/allWindowNoDup_flipRefPrunedSorted.fam
join -a 1 -a 2 -o0,2.2 -e ' -' $outputFolder/allWindowNoDup_flipRefPrunedSorted.fam $outputFolder/temp.pop > $outputFolder/temp2.pop
awk '{ print $2 }' $outputFolder/temp2.pop > $outputFolder/allWindowNoDup_flipRefPruned.pop


# /storage/software/admixture_linux-1.3.0/admixture $outputFolder/plink13/allWindowNoDup_flipRef4 --supervised
/storage/software/admixture_linux-1.3.0/admixture -s 314161 -j2 --supervised $outputFolder/allWindowNoDup_flipRefPruned.bed 2