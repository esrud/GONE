#script_genetic_distances.sh

########################################################

#Check number of arguments
if [ $# -ne 1 ]  
then
	echo "Usage: $0 <FILE>" 
	exit 1
fi

### Set arguments

FILE=$1  ### Data file name for files .ped and .map

###################### FILES NEEDED ########################

### data.ped
### data.map

### EXECUTABLES FILES NEEDED IN DIRECTORY PROGRAMMES:

### plink (needs gsl/2.1)

################### Obtain number of chromosomes ##################

NCHR=$(tail -n 1 $FILE.map | tr '\t' ' ' | cut -d ' ' -f1)

################### Divide file into chromosomes ##################
### The data.map and data.ped will be divided into the required chromosomes

for ((i=1; i<=$NCHR; i++))
do
./PROGRAMMES/plink --file $FILE --chr-set $NCHR --make-bed --chr $i --recode --out chromosome$i
done

################### INCLUDE GENETIC DISTANCES ##################
### Only for data for which genetic distances are available in the directory GENETIC_MAP,
### Distances between SNPs in cM are taken from directory GENETIC_MAP
### and included in .map file for each chromosome

for ((i=1; i<=$NCHR; i++))
do

./PROGRAMMES/plink --bfile chromosome$i --cm-map GENETIC_MAP/genetic_map_chr@.txt --make-bed --out chromosome_cms$i
./PROGRAMMES/plink --bfile chromosome_cms$i --recode --out chrom_cmsN$i

mv chrom_cmsN$i.map chromosome$i.map

paste chromosome$i.map >> KK3

done

for ((i=1; i<=$NCHR; i++))
do
rm chromosome_cms$i*
rm chrom_cmsN$i*
rm chromosome$i*
done

mv KK3 $FILE.map

###################################################################















