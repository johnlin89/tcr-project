# for all 129 samples
# for K in 1 2 3 4 5 6 7 8 9 10
# 	do 
# 		../admixture_macosx-1.3.0/admixture --cv allWindow_flip.bed $K | tee log${K}.out
# 	done
# 	
# grep -h CV log*.out
# K = 6 shows lowest error
outputFolder="/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults"

/storage/software/admixture_linux-1.3.0/admixture $outputFolder