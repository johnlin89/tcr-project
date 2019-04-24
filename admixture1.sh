# cross validation to determine right K
for K in 1 2 3 4 5
	do 
		../admixture_macosx-1.3.0/admixture --cv windowCovar_flip.bed $K | tee log${K}.out
	done