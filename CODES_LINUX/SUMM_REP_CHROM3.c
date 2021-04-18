// SUMM_REP_CHROM3.c

#include "libhdr"
#define BB 1001

int x, a, c, b, nchrom, NGEN, NBIN, PHASE, SAM, totSNP, numSNP;
unsigned long long int sumi[BB];
double w, sumcc[BB], sumd2[BB], sumD2[BB], sumVV[BB];
double cc[BB], d2[BB], D2[BB], VV[BB], F, G, numNIND;

FILE *fin, *fout, *fhwd, *fnsnp;

main()
{
	getintandskip("NGEN :",&NGEN, 1, 5000);
	getintandskip("NBIN :",&NBIN, 1, 1000);
	getintandskip("n.chrom :",&nchrom, 1, 301);

	/* ***************** readfile ******************** */

	fout = fopen ("outfileLD","w");
	fhwd = fopen ("outfileHWD","w");
	fnsnp = fopen ("nsnp","w");
	fin = fopen ("CHROM","r");

	for (c=1; c<=nchrom; c++)
	{
		fscanf(fin,"%d", &x);
		PHASE = x;
		fscanf(fin,"%d", &x);
		SAM = x;
		fscanf(fin,"%d", &x);
		numSNP = x;
		totSNP += numSNP;

		fscanf(fin,"%lf", &w);
		if (w != -99)
		{
			numNIND += w*numSNP;
		}

		fscanf(fin,"%lf", &w);
		if (w != -99)
		{
			F += w*numSNP;
		}

		for (b=1; b<=NBIN; b++)
		{
			fscanf(fin,"%d", &a);
			sumi[b] += a;
			fscanf(fin,"%lf", &w);
			sumcc[b] += w*a;
			fscanf(fin,"%lf", &w);
			sumd2[b] += w*a;
			fscanf(fin,"%d", &x);
			fscanf(fin,"%lf", &w);
			sumD2[b] += w*a;
			fscanf(fin,"%lf", &w);
			sumVV[b] += w*a;
		}
	}
	for (b=1; b<=NBIN; b++)
	{
		cc[b] = sumcc[b] / sumi[b];
		d2[b] = sumd2[b] / sumi[b];
		D2[b] = sumD2[b] / sumi[b];
		VV[b] = sumVV[b] / sumi[b];
	}

	fclose(fin);

	fprintf(fnsnp,"%d\n", totSNP);

	/* ***** outfileHWD ********************************************************* */

	fprintf(fhwd,"%f       Hardy-Weinberg deviation (sample)\n", F/totSNP);
	fprintf(fhwd,"%f       Hardy-Weinberg deviation (population)\n", (1.0 + ((F/totSNP)*(2.0*SAM - 1.0))) / (2.0*SAM - 1.0 + (F/totSNP)));

	/* ***** outfileLD ********************************************************* */

	fprintf(fout,"%d       Phase (0: pseudohaploids; 1: known phase; 2: unknown phase)\n", PHASE);
	fprintf(fout,"%f       sample size (individuals; corrected for zeroes)\n", numNIND/totSNP);
	fprintf(fout,"%f       Hardy-Weinberg deviation\n", F/totSNP);
	for (b=1; b<=NBIN; b++)
	if (sumi[b] != 0)
	{
		if (b <= 5)	fprintf(fout,"%lld %f %f %d\n", sumi[b], cc[b], D2[b]/VV[b], 2*b);
		else	fprintf(fout,"%lld %f %f %d\n", sumi[b], cc[b], D2[b]/VV[b], ((NGEN/NBIN)*(b-5))+10);
	}
}
