# GONE
Scripts and programs referred to in the paper "Recent demographic history inferred by high-resolution analysis of linkage disequilibrium" by Enrique Santiago, Irene Novo, Antonio F. Pardiñas, María Saura, Jinliang Wang and Armando Caballero. 
Molecular Biology and Evolution, 2020

Download the whole directory (Linux or MacOSX) as is and read "SCRIPT ESTIMATION Ne PROCEDURE.pdf"


# MODIFICATIONS
01/07/2020   Minor correction in LD_SNP_REAL3.c

07/07/2020   Minor modification in MANAGE_CHROMOSOMES2.c, INPUT_PARAMETERS_FILE and SCRIPT ESTIMATION Ne PROCEDURE.docx to note that a maximum of 100,000 SNPs are allowed per chromosome. A maximum number of 200 chromosomes and a maximum number of 1800 individuals are allowed.

25/08/2020   Minor modification in INPUT_PARAMETERS_FILE and SCRIPT ESTIMATION Ne PROCEDURE.docx to note that the option cMMb=0 (for human data) does not hold. If the .map file has genetic distances (third column), these will be used. If that column has zeroes (no genetic distances are available) then an average rate of recombination of cMMb (as included in the INPUT_PARAMETERS_FILE) will be assumed.

07/09/2020   Minor modification in LD_SNP_REAL3.c so that the string "-9" in the phenotypic column of the ped file is not necessary anymore and that string can be present in the pedfile without causing trouble.

13/09/2020   Corrected an error in the function to get random numbers in MANAGE_CHROMOSOMES2.c. Now the file called seedfile is necessary in the running directory to start the random sampling of SNPs within each chromosome. Every time a run is finished the seedfile will be changed automatically.
