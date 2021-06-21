# GONE
Scripts and programs referred to in the paper "Recent demographic history inferred by high-resolution analysis of linkage disequilibrium" by Enrique Santiago, Irene Novo, Antonio F. Pardiñas, María Saura, Jinliang Wang and Armando Caballero. Molecular Biology and Evolution, 2020 Volume 37, Issue 12, Pages 3642–3653, https://doi.org/10.1093/molbev/msaa169

# USE
To use the software please download the whole appropriate directory (Linux or MacOSX). These directories include the scripts and executables files of all necessary programmes. See USER´S GUIDE.
The codes of the programmes are available in the directories CODES, but they are not necessary to run the software. 

# MODIFICATIONS
01/07/2020   Minor correction in LD_SNP_REAL3.c

07/07/2020   Minor modification in MANAGE_CHROMOSOMES2.c, INPUT_PARAMETERS_FILE and SCRIPT ESTIMATION Ne PROCEDURE.docx to note that a maximum of 100,000 SNPs are allowed per chromosome. A maximum number of 200 chromosomes and a maximum number of 1800 individuals are allowed.

25/08/2020   Minor modification in INPUT_PARAMETERS_FILE and SCRIPT ESTIMATION Ne PROCEDURE.docx to note that the option cMMb=0 (for human data) does not hold. If the .map file has genetic distances (third column), these will be used. If that column has zeroes (no genetic distances are available) then an average rate of recombination of cMMb (as included in the INPUT_PARAMETERS_FILE) will be assumed.

07/09/2020   Minor modification in LD_SNP_REAL3.c so that the string "-9" in the phenotypic column of the ped file is not necessary anymore and that string can be present in the pedfile without causing trouble.

13/09/2020   Corrected an error in the function to get random numbers in MANAGE_CHROMOSOMES2.c. Now the file called seedfile is necessary in the running directory to start the random sampling of SNPs within each chromosome. Every time a run is finished the seedfile will be changed automatically.

15/10/2020   Small correction in SCRIPT ESTIMATION Ne PROCEDURE.docx

26/10/2020   Small correction in script_GONE.sh to generate random seed, and the corresponding modification in SCRIPT ESTIMATION Ne PROCEDURE.pdf

09/01/2021   Small addition in SCRIPT ESTIMATION Ne PROCEDURE.docx: NOTE: If the population has recent migrants from another population, the estimation of Ne will be biased. A typical artefact observed is a very recent drastic drop and a previous increase (see Fig. 2f of manuscript). This can be partly corrected by using a maximum value of c lower than that recommended above, for example hc=0.01. 

18/04/2021   Modification of programmes LD_SNP_REAL3.c and SUMM_REP_CHROM3.c to calculate the deviations from Hardy-Weinberg proportions (Wrigth's Fis). Estimates for the sample and for the population are now shown in the outfileHWD output file. If the estimate for the population substantialy deviates from the expected 0 value of a panmictic population (say Fis > 0.02 or Fis < -0.02), this may imply certain adxmixture in your sample and some artefacts, tipically a sudden recent drop in Ne, may be found.

21/06/2021   New USER´S GUIDE which substitutes the old SCRIPT ESTIMATION Ne PROCEDURE.docx
