outputFolder="/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink15"
ref_dir="/storage/mips/MIPS_Updated.2019-02-21/data/1KG_Phase3"

# supervised learning

# # create .pop file
# awk '{ if (NR!=1) {print $1 " " $6} }' $ref_dir/CEU-YRI.psam > $outputFolder/temp.pop
# sort $outputFolder/temp.pop -o $outputFolder/temp.pop
# sort $outputFolder/allWindowNoDup_flipRefPruned.fam -o $outputFolder/allWindowNoDup_flipRefPrunedSorted.fam
# join -a 1 -a 2 -o0,2.2 -e ' -' $outputFolder/allWindowNoDup_flipRefPrunedSorted.fam $outputFolder/temp.pop > $outputFolder/temp2.pop
# awk '{ print $2 }' $outputFolder/temp2.pop > $outputFolder/allWindowNoDup_flipRefPruned.pop


# # /storage/software/admixture_linux-1.3.0/admixture $outputFolder/plink13/allWindowNoDup_flipRef4 --supervised
# /storage/software/admixture_linux-1.3.0/admixture -B -s 314161 -j2 --supervised $outputFolder/allWindowNoDup_flipRefPruned.bed 2

# mv allWindowNoDup_flipRefPruned.2.* ../admixtureResults/admixture1/

# unsupervised learning for all 336 samples
# for K in 1 2 3 4 5
# 	do 
# 		/storage/software/admixture_linux-1.3.0/admixture -s 314161 -j8 --cv=10 $outputFolder/allWindowNoDup_flipRefPruned.bed  $K | tee log${K}.out
# 	done

# grep -h CV log*.out > crossValidation.txt
# mv allWindowNoDup_flipRefPruned* ../admixtureResults/admixture2/
# mv crossValidation.txt ../admixtureResults/admixture2/

# unsupervised learning for all 129 samples
# for K in 1 2 3 4 5
# 	do 
# 		/storage/software/admixture_linux-1.3.0/admixture -s 314161 -j8 --cv=10 /storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink8/allWindow.bed  $K | tee log${K}.out
# 	done

# grep -h CV log*.out > crossValidation.txt
# mv allWindow* ../admixtureResults/admixture3/
# mv crossValidation.txt ../admixtureResults/admixture3/

/storage/software/admixture_linux-1.3.0/admixture -B -s 314161 -j8 /storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink8/allWindow.bed 2

# K = 6 shows lowest error