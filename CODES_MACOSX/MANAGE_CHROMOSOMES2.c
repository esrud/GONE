
// MANAGE_CHROMOSOMES2.c

#include "libhdr"
#define CC 200
#define SS 1000000

int maxNSNP, x, i, c, s, k, NCHR, NIND, snp_nchrom[CC], ranSNP[CC][SS];
double w;
char S[CC], ch;

FILE *fnchr, *fnind, *fsnpchrom, *fmap, *fped, *fcheck;

main()
{
//	fcheck = fopen ("checkfile","w");

	getseed();
	getintandskip("maxNSNP :",&maxNSNP, -99, 100000);

	// ******************** Read NCHR, NIND and NSNP ******************** 

	fnchr = fopen ("NCHR","r");
	fscanf(fnchr,"%d", &x);
	NCHR = x;
	// printf("NCHR=%d\n\n", NCHR);

	fnind = fopen ("NIND","r");
	fscanf(fnind,"%d", &x);
	NIND = x;
	// printf("NIND=%d\n\n", NIND);

	fsnpchrom = fopen ("SNP_CHROM","r");
	for (c=1; c<=NCHR; c++)
	{
		fscanf(fsnpchrom,"%d", &x);
		snp_nchrom[c] = x;
		// if (c==1) printf("NSNP_CHROM[%d]=%d\n", c, snp_nchrom[c]); 
	}

	// ******** Randomize SNPs *******

	if (maxNSNP != -99)
	for (c=1; c<=NCHR; c++)
	{
		if (maxNSNP >= snp_nchrom[c])
		{
			for (s=1; s<=snp_nchrom[c]; s++)
			ranSNP[c][s]=1;
		}
		else
		{
			for (s=1; s<=snp_nchrom[c]; s++)
			if (uniform() <= ((double)maxNSNP/(double)snp_nchrom[c])) ranSNP[c][s]=1;
		}
	}

	// ******************** Read data.map ******************** 

	FILE *mapchrom[NCHR];
	for (c=1; c<=NCHR; c++)
	{
		char filename[20];
		sprintf(filename, "chromosome%d.map", c);
		mapchrom[c] = fopen(filename, "w");
	}

	fmap = fopen ("data.map","r");

	for (c=1; c<=NCHR; c++)
	{
		for (s=1; s<=snp_nchrom[c]; s++)
		{
			fscanf(fmap,"%d", &x);
			if (maxNSNP == -99)	fprintf(mapchrom[c], "%d\t", x);
			else			if (ranSNP[c][s] == 1)	fprintf(mapchrom[c], "%d\t", x);

			fscanf(fmap,"%s", &S);
			if (maxNSNP == -99)	fprintf(mapchrom[c], "%s\t", S);
			else			if (ranSNP[c][s] == 1)	fprintf(mapchrom[c], "%s\t", S);

			fscanf(fmap,"%lf", &w);
			if (maxNSNP == -99)	fprintf(mapchrom[c], "%f\t", w);
			else			if (ranSNP[c][s] == 1)	fprintf(mapchrom[c], "%f\t", w);

			fscanf(fmap,"%d", &x);
			if (maxNSNP == -99)	fprintf(mapchrom[c], "%d\n", x);
			else			if (ranSNP[c][s] == 1)	fprintf(mapchrom[c], "%d\n", x);
		}
	}

	fclose(fmap);

	// ******************** Read data.ped ******************** 

	FILE *pedchrom[NCHR];
	for (c=1; c<=NCHR; c++)
	{
		char filename[20];
		sprintf(filename, "chromosome%d.ped", c);
		pedchrom[c] = fopen(filename, "w");
	}

	fped = fopen ("data.ped","r");

	for (i=1; i<=NIND; i++)
	{
		fscanf(fped,"%s", &S);
	//	if ((i==1)||(i==2)) fprintf(fcheck, "%s ", S);
		for (c=1; c<=NCHR; c++)	fprintf(pedchrom[c], "%s ", S);

		fscanf(fped,"%s", &S);
	//	if ((i==1)||(i==2)) fprintf(fcheck, "%s ", S);
		for (c=1; c<=NCHR; c++)	fprintf(pedchrom[c], "%s ", S);

		fscanf(fped,"%s", &S);
	//	if ((i==1)||(i==2)) fprintf(fcheck, "%s ", S);
		for (c=1; c<=NCHR; c++)	fprintf(pedchrom[c], "%s ", S);

		fscanf(fped,"%s", &S);
	//	if ((i==1)||(i==2)) fprintf(fcheck, "%s ", S);
		for (c=1; c<=NCHR; c++)	fprintf(pedchrom[c], "%s ", S);

		fscanf(fped,"%s", &S);
	//	if ((i==1)||(i==2)) fprintf(fcheck, "%s ", S);
		for (c=1; c<=NCHR; c++)	fprintf(pedchrom[c], "%s ", S);

		fscanf(fped,"%s", &S);
	//	if ((i==1)||(i==2)) fprintf(fcheck, "-9 ");
		for (c=1; c<=NCHR; c++)	fprintf(pedchrom[c], "-9 "); 

		for (c=1; c<=NCHR; c++)
		{
			for (s=1; s<=snp_nchrom[c]; s++)
			{
				fscanf(fped,"%s", &S);
				if (maxNSNP == -99)	fprintf(pedchrom[c], "%s ", S);
				else			if (ranSNP[c][s] == 1)	fprintf(pedchrom[c], "%s ", S);			
			//	if ((i==1)||(i==2)) fprintf(fcheck, "%s ", S);

				fscanf(fped,"%s", &S);
				if (maxNSNP == -99)	fprintf(pedchrom[c], "%s ", S);
				else			if (ranSNP[c][s] == 1)	fprintf(pedchrom[c], "%s ", S);			
			//	if ((i==1)||(i==2)) fprintf(fcheck, "%s ",S);
			}
			fprintf(pedchrom[c], "\n");
		}
	}

	fclose(fped);
	writeseed();
	
	return(0);
}

/* **************************************************************************** */
