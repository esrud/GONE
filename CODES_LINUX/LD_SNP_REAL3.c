
// LD_SNP_REAL3.c (LD estimates from data.ped, with and without PHASE)

#include "libhdr"

#define II 1801 // Maximum N=1800 individuals
#define SS 120001 // Maximum 120000 SNPs segregating
#define BB 1001   // Maximum number of bins (1000)

unsigned long long int n_[BB];
int n, b, NGEN, NBIN, snp_m2, snp_mZ, snp_mono, ZERO, maxCEROS, s1, s2, DISTANCE;
int i, j, s, ss, k, x, NIND, NSNP, snp[II][SS][2], all[SS][5], pos[SS];
int PHASE, SUM_ALLELES[SS], numSNP, numNIND, calc_SNP[SS], mono_SNP[SS], allele[SS];
int gam11, gam12, gam21, gam22;

double c, g, MAF, mor, dis[SS], cM, TOTAL_DIS;
double Ne_r2_[BB], D2_[BB], VV_[BB], d2_[BB], d2_s_[BB], Ne_d2_[BB];
double freq[SS][3], freqs, freqss, freqG11s, freqG11ss, freqG1111, freqG1112, freqG1121, freqG1212;
double D, D2, VV, r2, w, r2_[BB], cc, cc_[BB], qi, freq_p2s, freq_MAF, f, fp;

char ch, S[100], filename[20];

char *x1, *x2, *x3, *x4, *x5, *x6, *x7, *x8, *x9;

struct acc AVE_D_[BB], AVE_D2_[BB], AVE_VV_[BB], AVE_r2_[BB], AVE_c_[BB], AVE_numNIND;
struct acc P_p2, pq;

FILE *fped, *fmap, *fout, *fpar;

int main(int argc, char *argv[])
{
	s1 = 1;
	s2 = 5;
	tracelevel=0;

	x1 = argv[1];
	x2 = argv[2];
	x3 = argv[3];
	x4 = argv[4];
	x5 = argv[5];
	x6 = argv[6];
	x7 = argv[7];
	x8 = argv[8];
	x9 = argv[9];

	n = atoi(x1);
	NIND = atoi(x2);
	MAF = atof(x3);
	PHASE = atoi(x4);
	NGEN = atoi(x5);
	NBIN = atoi(x6);
	ZERO = atoi(x7);

	if (ZERO==0)  maxCEROS = 0;
	else 	    maxCEROS = 1000;

	DISTANCE = atoi(x8);
	cM = atof(x9);

	sprintf(filename, "outfileLD%d", n);
	fout = fopen (filename,"w");

	sprintf(filename, "parameters%d", n);
	fpar = fopen (filename,"w");

	readfiles();
	LD();

	return(0);
}

/* **************************************************************************** */

readfiles()
{
	// ********** Read data.map file to get SNP positions ********** 

	sprintf(filename, "chromosome%d.map", n);
	fmap = fopen (filename,"r");

	while (!feof(fmap))
	{
		NSNP ++;
		fscanf(fmap,"%d", &x);
		fscanf(fmap,"%s", &S);

		fscanf(fmap,"%lf", &w);
		dis[NSNP] = w;
		TOTAL_DIS += dis[NSNP];
		fscanf(fmap,"%d", &x);
		pos[NSNP] = x;
//		if (tracelevel != 0)	fprintf(fout, "dis=%f  pos=%d  ", dis[NSNP], pos[NSNP]);
	}
	NSNP = NSNP-1;
//	if (tracelevel != 0)	fprintf(fout, "\n\nNSNP=%d\n\n", NSNP);

	fclose(fmap);

	// ********** Read data.ped file to get genotypes ********** 

	sprintf(filename, "chromosome%d.ped", n);
	fped = fopen (filename,"r");

	for (i=1; i<=NIND; i++)
	{
		for (j=1; j<=6; j++)
		fscanf(fped,"%s", &S);

		for (s=1; s<=NSNP; s++)
		{
			fscanf(fped,"%c", &ch);
			fscanf(fped,"%c", &ch);
			if (ch == '1')	ch = 'A';
			if (ch == '2')	ch = 'T';

			if (ch == 'A')	snp[i][s][1] = 1;
			else if (ch == 'T')	snp[i][s][1] = 2;
			else if (ch == 'C')	snp[i][s][1] = 3;
			else if (ch == 'G')	snp[i][s][1] = 4;
			else if (ch == '0')	snp[i][s][1] = 0;

			fscanf(fped,"%c", &ch);
			fscanf(fped,"%c", &ch);
			if (ch == '1')	ch = 'A';
			if (ch == '2')	ch = 'T';

			if (ch == 'A')	snp[i][s][2] = 1;
			else if (ch == 'T')	snp[i][s][2] = 2;
			else if (ch == 'C')	snp[i][s][2] = 3;
			else if (ch == 'G')	snp[i][s][2] = 4;
			else if (ch == '0')	snp[i][s][2] = 0;
		}
	}

	if (tracelevel != 0)
	{
		for (s=1; s<=NSNP; s++)
		{
			if (s == s1)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\n");
			}
			if (s == s2)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\n");
			}
		}
	}

	fclose(fped);

	return(0);
}

/* **************************************************************************** */

LD()
{
	/******* Alleles segregating *******/

	for (s=1; s<=NSNP; s++)
	for (i=1; i<=NIND; i++)
	for (k=1; k<=2; k++)
	{
		if      (snp[i][s][k] == 1) 	all[s][1] ++;
		else if (snp[i][s][k] == 2) 	all[s][2] ++;
		else if (snp[i][s][k] == 3) 	all[s][3] ++;
		else if (snp[i][s][k] == 4) 	all[s][4] ++;
		else if (snp[i][s][k] == 0) 	all[s][5] ++;
	}
	for (s=1; s<=NSNP; s++)	SUM_ALLELES[s] = (all[s][1] + all[s][2] + all[s][3] + all[s][4] + all[s][5]);

	for (s=1; s<=NSNP; s++)
	{
		if (all[s][1] > 0) 	allele[s]++;
		if (all[s][2] > 0) 	allele[s]++;
		if (all[s][3] > 0) 	allele[s]++;
		if (all[s][4] > 0) 	allele[s]++;
	}

	if (tracelevel != 0)		
	{
		fprintf(fout, "\n\n");
		for (s=1; s<=NSNP; s++)	if ((s == s1)||(s == s2)) fprintf(fout, "s=%d  %d %d %d %d %d\n", s, all[s][1], all[s][2], all[s][3], all[s][4], all[s][5]);
//		for (s=1; s<=NSNP; s++)	if (allele[s]>2) fprintf(fout, "s=%d  %d %d %d %d %d\n", s, all[s][1], all[s][2], all[s][3], all[s][4], all[s][5]);
		fprintf(fout, "\n\n");
	}

	/******* Change allele codes to 1 and 2 *******/

	for (s=1; s<=NSNP; s++)
	if (allele[s] <= 2)
	if (all[s][5] <= maxCEROS)
	{
		if ( (all[s][1] != 0) && (all[s][2] != 0) )
		{
			freq[s][1] = all[s][1]/((2.0*NIND)-all[s][5]);
			freq[s][2] = all[s][2]/((2.0*NIND)-all[s][5]);
			calc_SNP[s] = 1;
			numSNP++;
		}
		else if ( (all[s][1] != 0) && (all[s][3] != 0) )
		{
			freq[s][1] = all[s][1]/((2.0*NIND)-all[s][5]);
			freq[s][2] = all[s][3]/((2.0*NIND)-all[s][5]);

			for (i=1; i<=NIND; i++)
			for (k=1; k<=2; k++)
			if (snp[i][s][k] == 3) 	snp[i][s][k] = 2;
			calc_SNP[s] = 1;
			numSNP++;
		}
		else if ( (all[s][1] != 0) && (all[s][4] != 0) )
		{
			freq[s][1] = all[s][1]/((2.0*NIND)-all[s][5]);
			freq[s][2] = all[s][4]/((2.0*NIND)-all[s][5]);

			for (i=1; i<=NIND; i++)
			for (k=1; k<=2; k++)
			if (snp[i][s][k] == 4) 	snp[i][s][k] = 2;
			calc_SNP[s] = 1;
			numSNP++;
		}
		else if ( (all[s][2] != 0) && (all[s][3] != 0) )
		{
			freq[s][1] = all[s][2]/((2.0*NIND)-all[s][5]);
			freq[s][2] = all[s][3]/((2.0*NIND)-all[s][5]);

			for (i=1; i<=NIND; i++)
			for (k=1; k<=2; k++)
			{
				if (snp[i][s][k] == 2) 	snp[i][s][k] = 1;
				if (snp[i][s][k] == 3) 	snp[i][s][k] = 2;
			}
			calc_SNP[s] = 1;
			numSNP++;
		}
		else if ( (all[s][2] != 0) && (all[s][4] != 0) )
		{
			freq[s][1] = all[s][2]/((2.0*NIND)-all[s][5]);
			freq[s][2] = all[s][4]/((2.0*NIND)-all[s][5]);

			for (i=1; i<=NIND; i++)
			for (k=1; k<=2; k++)
			{
				if (snp[i][s][k] == 2) 	snp[i][s][k] = 1;
				if (snp[i][s][k] == 4) 	snp[i][s][k] = 2;
			}
			calc_SNP[s] = 1;
			numSNP++;
		}
		else if ( (all[s][3] != 0) && (all[s][4] != 0) )
		{
			freq[s][1] = all[s][3]/((2.0*NIND)-all[s][5]);
			freq[s][2] = all[s][4]/((2.0*NIND)-all[s][5]);

			for (i=1; i<=NIND; i++)
			for (k=1; k<=2; k++)
			{
				if (snp[i][s][k] == 3) 	snp[i][s][k] = 1;
				if (snp[i][s][k] == 4) 	snp[i][s][k] = 2;
			}
			calc_SNP[s] = 1;
			numSNP++;
		}
		else
		{
			snp_mono ++;
			mono_SNP[s] = 1;
//			printf("s=%d  ", s);
		}
	}

	if (tracelevel != 0)	
	{
/*		fprintf(fout, "\MONO, MORE THAN 2 alleles AND MORE THAN Maxceros\n");
		for (s=1; s<=NSNP; s++)
		{
			if (mono_SNP[s] == 1)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
			if (allele[s] > 2)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
			if (all[s][5] > maxCEROS)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
		}
		fprintf(fout, "\nCHANGED\n");
*/
		for (s=1; s<=NSNP; s++)
		{
			if (s == s1)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
			if (s == s2)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
		}
	}

	/******* Major allele is 1 and minor is 2 *******/

	for (s=1; s<=NSNP; s++)
	if (calc_SNP[s] == 1)
	{
		if (freq[s][1] < freq[s][2])
		{
			for (i=1; i<=NIND; i++)
			for (k=1; k<=2; k++)
			{
				if (snp[i][s][k] == 1) 	snp[i][s][k] = 2;
				else if (snp[i][s][k] == 2) snp[i][s][k] = 1;
			}
			w=freq[s][1]; freq[s][1]=freq[s][2]; freq[s][2]=w;
		}
	}

	if (tracelevel != 0)	
	{
		fprintf(fout, "\nCHANGED major is 1\n");
		for (s=1; s<=NSNP; s++)
		{
			if (s == s1)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
			if (s == s2)
			{
				fprintf(fout, "s=%d  ", s);
				for (i=1; i<=NIND; i++)	 fprintf(fout, "%d%d ", snp[i][s][1], snp[i][s][2]);
				fprintf(fout, "\nqs=%f ps=%f\n\n", freq[s][1], freq[s][2]);
			}
		}
	}

	/******* Calculation of D, D2, VV and r2 *******/

	freq_MAF=1.0;

	for (s=1; s<NSNP; s++)
	for (ss=s+1; ss<=NSNP; ss++)
	if ((allele[s] <= 2) && (allele[ss] <= 2))
	if ((calc_SNP[s] == 1) && (calc_SNP[ss] == 1))
	if ((freq[s][2] > MAF) && (freq[ss][2] > MAF))
	{
		freqs=0.0; freqss=0.0; freqG11s=0.0; freqG11ss=0.0;
		freqG1111=0.0; freqG1112=0.0; freqG1121=0.0; freqG1212=0.0;
		freq_p2s=0.0;

		gam11=0; gam12=0; gam21=0; gam22=0; 

		if ((ss == s+1) && (freq[s][2] < freq_MAF))	freq_MAF = freq[s][2];

		numNIND=0;
		for (i=1; i<=NIND; i++)
		if ( (snp[i][s][1] != 0) && (snp[i][s][2] != 0) && (snp[i][ss][1] != 0) & (snp[i][ss][2] != 0) )
		{
			numNIND ++;

			if ((snp[i][s][1] == 1) && (snp[i][s][2] == 1))	qi = 1.0;
			if ((snp[i][s][1] == 1) && (snp[i][s][2] == 2))	qi = 0.5;
			if ((snp[i][s][1] == 2) && (snp[i][s][2] == 1))	qi = 0.5;
			if ((snp[i][s][1] == 2) && (snp[i][s][2] == 2))	qi = 0.0;

			freq_p2s += (qi * qi);

			if (snp[i][s][1] == 1) freqs ++;
			if (snp[i][s][2] == 1) freqs ++;

			if (snp[i][ss][1] == 1) freqss ++;
			if (snp[i][ss][2] == 1) freqss ++;

			if ( (snp[i][s][1] == 1) && (snp[i][s][2] == 1) )
		 	freqG11s ++;

			if ( (snp[i][ss][1] == 1) && (snp[i][ss][2] == 1) )
		 	freqG11ss ++;

			if ( (snp[i][s][1] == 1) && (snp[i][ss][1] == 1) && (snp[i][s][2] == 1) && (snp[i][ss][2] == 1) )
		 	freqG1111 ++;

			if ( (snp[i][s][1] == 1) && (snp[i][ss][1] == 1) && (snp[i][s][2] == 1) && (snp[i][ss][2] == 2) )
		 	freqG1112 ++;
			if ( (snp[i][s][1] == 1) && (snp[i][ss][1] == 2) && (snp[i][s][2] == 1) && (snp[i][ss][2] == 1) )
		 	freqG1112 ++;

			if ( (snp[i][s][1] == 1) && (snp[i][ss][1] == 1) && (snp[i][s][2] == 2) && (snp[i][ss][2] == 1) )
		 	freqG1121 ++;
			if ( (snp[i][s][1] == 2) && (snp[i][ss][1] == 1) && (snp[i][s][2] == 1) && (snp[i][ss][2] == 1) )
		 	freqG1121 ++;

			if ( (snp[i][s][1] == 1) && (snp[i][ss][1] == 1) && (snp[i][s][2] == 2) && (snp[i][ss][2] == 2) )
		 	freqG1212 ++;
			if ( (snp[i][s][1] == 2) && (snp[i][ss][1] == 2) && (snp[i][s][2] == 1) && (snp[i][ss][2] == 1) )
		 	freqG1212 ++;
			if ( (snp[i][s][1] == 1) && (snp[i][ss][1] == 2) && (snp[i][s][2] == 2) && (snp[i][ss][2] == 1) )
		 	freqG1212 ++;
			if ( (snp[i][s][1] == 2) && (snp[i][ss][1] == 1) && (snp[i][s][2] == 1) && (snp[i][ss][2] == 2) )
		 	freqG1212 ++;

			// PHASE 1
			if      ( (snp[i][s][1] == 1) && (snp[i][ss][1] == 1) ) 	gam11 ++;
			else if ( (snp[i][s][1] == 1) && (snp[i][ss][1] == 2) ) 	gam12 ++;
			else if ( (snp[i][s][1] == 2) && (snp[i][ss][1] == 1) ) 	gam21 ++;
			else if ( (snp[i][s][1] == 2) && (snp[i][ss][1] == 2) ) 	gam22 ++;

			if      ( (snp[i][s][2] == 1) && (snp[i][ss][2] == 1) ) 	gam11 ++;
			else if ( (snp[i][s][2] == 1) && (snp[i][ss][2] == 2) ) 	gam12 ++;
			else if ( (snp[i][s][2] == 2) && (snp[i][ss][2] == 1) ) 	gam21 ++;
			else if ( (snp[i][s][2] == 2) && (snp[i][ss][2] == 2) ) 	gam22 ++;
		}

	       if (numNIND > NIND/2)
	       {

		accum(&AVE_numNIND, 1.0/numNIND);

		freq_p2s = freq_p2s/numNIND;
		freqs = freqs/(2.0*numNIND);
		freqss = freqss/(2.0*numNIND);
		freqG11s = freqG11s/numNIND;
		freqG11ss = freqG11ss/numNIND;
		freqG1111 = freqG1111/numNIND;
		freqG1112 = freqG1112/numNIND;
		freqG1121 = freqG1121/numNIND;
		freqG1212 = freqG1212/numNIND;

		if (tracelevel != 0)	
		{
			if ((s == s1) && (ss == s2))
			{
				fprintf(fout, "\ns=%d  ss=%d  numNIND=%d\n", s, ss, numNIND);
				fprintf(fout, "freqs=%f  freqss=%f\n", freqs, freqss);
				fprintf(fout, "freqG11s=%f  freqG11ss=%f\n", freqG11s, freqG11ss);
				fprintf(fout, "freqG1111=%f\n", freqG1111);
				fprintf(fout, "freqG1112=%f\n", freqG1112);
				fprintf(fout, "freqG1121=%f\n", freqG1121);
				fprintf(fout, "freqG1212=%f\n", freqG1212);
				fprintf(fout, "gam11=%d gam12=%d gam21=%d gam22=%d\n", gam11, gam12, gam21, gam22);
				fprintf(fout, "\n");
			}
/*			if (numNIND < NIND)
			{
				fprintf(fout, "\ns=%d  ss=%d  numNIND=%d\n", s, ss, numNIND);
				fprintf(fout, "freqs=%f  freqss=%f\n", freqs, freqss);
				fprintf(fout, "freqG11s=%f  freqG11ss=%f\n", freqG11s, freqG11ss);
				fprintf(fout, "freqG1111=%f\n", freqG1111);
				fprintf(fout, "freqG1112=%f\n", freqG1112);
				fprintf(fout, "freqG1121=%f\n", freqG1121);
				fprintf(fout, "freqG1212=%f\n", freqG1212);
				fprintf(fout, "\n");
			}
*/		}

		if ((PHASE == 0)||(PHASE == 1))		D = ( (gam11/(2.0*numNIND)) * (gam22/(2.0*numNIND)) ) - ( (gam12/(2.0*numNIND)) * (gam21/(2.0*numNIND)) );
		else if (PHASE == 2)	D = (2.0*freqG1111) + freqG1112 + freqG1121 + (freqG1212/2.0) - (2.0*freqs*freqss);

		D2 = (D * D);
		VV = ( freqs * (1.0 - freqs) ) * ( freqss * (1.0 - freqss) );
		r2 = D2 / VV;


		if (TOTAL_DIS == 0.0)	mor = ( abs(pos[ss]-pos[s]) / (1000000/cM) ) * 0.01;
		else					mor = (dis[ss]-dis[s])*0.01;

		if (DISTANCE == 0)	c = mor;
		else	if (DISTANCE == 1)	c = (1.0 - exp(-2.0*mor))/2.0;
		else	if (DISTANCE == 2)	c = 0.5*(exp(4.0*mor)-1)/(exp(4.0*mor)+1);
		if (c>0.5) c=0.5;
		g = 1.0/(2.0*c);

		if (tracelevel != 0)	
		{
			if ((s == s1) && (ss == s2)) fprintf(fout, "\n\nSNPs %d-%d  qs=%f ps=%f qss=%f pss=%f D=%f D2=%f VV=%f r2=%f c=%f g=%f\n\n",
				s, ss, freqs, 1.0-freqs, freqss, 1.0-freqss, D, D2, VV, r2, c, g);
		}

		if (tracelevel != 0)
		{
			fprintf(fout, "SNPs %d-%d  qs=%f ps=%f qss=%f pss=%f D=%f D2=%f VV=%f r2=%f c=%f g=%f\n",
				s, ss, freqs, 1.0-freqs, freqss, 1.0-freqss, D, D2, VV, r2, c, g);
		}

			if ( (freqs*(1.0-freqs) != 0.0) && (freqss*(1.0-freqss) != 0.0) )
			{
				for (b=1; b<=NBIN; b++)
				if (c != 0.0)
				{
					if (b==1)
					{
						if ((g > 0) && (g <= 2))
						{
							n_[b] ++;
							accum(&AVE_D_[b], D);
							accum(&AVE_D2_[b], D2);
							accum(&AVE_VV_[b], VV);
							accum(&AVE_r2_[b], r2);
							accum(&AVE_c_[b], 1.0/c);
						}
					}
					else if (b==2)
					{
						if ((g > 2) && (g <= 4))
						{
							n_[b] ++;
							accum(&AVE_D_[b], D);
							accum(&AVE_D2_[b], D2);
							accum(&AVE_VV_[b], VV);
							accum(&AVE_r2_[b], r2);
							accum(&AVE_c_[b], 1.0/c);
						}
					}
					else if (b==3)
					{
						if ((g > 4) && (g <= 6))
						{
							n_[b] ++;
							accum(&AVE_D_[b], D);
							accum(&AVE_D2_[b], D2);
							accum(&AVE_VV_[b], VV);
							accum(&AVE_r2_[b], r2);
							accum(&AVE_c_[b], 1.0/c);
						}
					}
					else if (b==4)
					{
						if ((g > 6) && (g <= 8))
						{
							n_[b] ++;
							accum(&AVE_D_[b], D);
							accum(&AVE_D2_[b], D2);
							accum(&AVE_VV_[b], VV);
							accum(&AVE_r2_[b], r2);
							accum(&AVE_c_[b], 1.0/c);
						}
					}
					else if (b==5)
					{
						if ((g > 8) && (g <= 10))
						{
							n_[b] ++;
							accum(&AVE_D_[b], D);
							accum(&AVE_D2_[b], D2);
							accum(&AVE_VV_[b], VV);
							accum(&AVE_r2_[b], r2);
							accum(&AVE_c_[b], 1.0/c);
						}
					}
					else
					{
						if ( (g > ((NGEN/NBIN)*(b-6))+10) && (g <= ((NGEN/NBIN)*(b-5))+10) )
						{
							n_[b] ++;
							accum(&AVE_D_[b], D);
							accum(&AVE_D2_[b], D2);
							accum(&AVE_VV_[b], VV);
							accum(&AVE_r2_[b], r2);
							accum(&AVE_c_[b], 1.0/c);
						}
					}
				}
			}
	       }
	}

	/******* Calculation of Hardy-Weinberg deviation *******/

	for (s=1; s<=NSNP; s++)
	if (calc_SNP[s] == 1)
	{
		freqs=0.0; freqG11s=0.0;

		numNIND=0;
		for (i=1; i<=NIND; i++)
		if ((snp[i][s][1] != 0) && (snp[i][s][2] != 0))
		{
			numNIND ++;

			if (snp[i][s][1] == 1) freqs ++;
			if (snp[i][s][2] == 1) freqs ++;

			if ( (snp[i][s][1] == 1) && (snp[i][s][2] == 1) )
		 	freqG11s ++;
		}

		freqs = freqs/(2.0*numNIND);
		freqG11s = freqG11s/numNIND;

		if (freqs*(1.0-freqs) != 0.0)
		{
			accum(&P_p2, freqG11s-(freqs*freqs));
			accum(&pq, freqs*(1.0-freqs));

			if (tracelevel != 0)	
			{
				fprintf(fout, "\n****** s=%d *******\n", s);
				fprintf(fout, "freqs=%f  freqs*(1.0-freqs)=%f\n", freqs, freqs*(1.0-freqs));
				fprintf(fout, "freqG11s=%f  freqG11s-(freqs*freqs)=%f\n", freqG11s, freqG11s-(freqs*freqs));
				fprintf(fout, "\n");
			}
		}
	}

	/* **************************************** */

	// Table of results

	for (b=1; b<=NBIN; b++)
	{
		cc_[b] = 1.0/accmean(&AVE_c_[b]);
		r2_[b] = accmean(&AVE_r2_[b]);
		D2_[b] = accmean(&AVE_D2_[b]);
		VV_[b] = accmean(&AVE_VV_[b]);

//		if (NIND == 2)	d2_[b] = (D2_[b] / VV_[b]);
//		else		d2_[b] = ( (D2_[b] / (4.0*VV_[b])) );

		d2_[b] = (D2_[b] / VV_[b]);

		Ne_r2_[b] = (1.0/(4.0*cc_[b]*r2_[b])) - (1.0/(2.0*cc_[b]));
		Ne_d2_[b] = (1.0/(4.0*cc_[b]*d2_[b])) - (1.0/(2.0*cc_[b]));
	}

	for (s=1; s<=NSNP; s++)
	{
		if (allele[s] > 2)	snp_m2 ++;
		if (all[s][5] > maxCEROS)	snp_mZ ++;
	}

	f = accmean(&P_p2) / accmean(&pq);
	fp = (1.0 + (f*(2.0*NIND - 1.0))) / (2.0*NIND - 1.0 + f);

	fprintf(fpar," NIND(real sample)=%d\n NSNP=%d\n NSNP_calculations=%d\n NSNP_+2alleles=%d\n NSNP_zeroes=%d\n NSNP_monomorphic=%d\n NIND_corrected=%f\n freq_MAF=%f\n F_dev_HW (sample)=%f\n F_dev_HW (pop)=%f\n",
			NIND, NSNP, numSNP, snp_m2, snp_mZ, snp_mono, 1.0/accmean(&AVE_numNIND), freq_MAF, f, fp);

	if (TOTAL_DIS == 0.0)	fprintf(fpar," No genetic distances; using %f cM per Mb\n\n", cM);
	else					fprintf(fpar," Genetic distances available in map file\n\n");

//	OUTFILELD
//	if (tracelevel != 0)	fprintf(fout,"\nn   c    d2    g    D2   VV\n");

	fprintf(fout,"%d\n", PHASE);
	fprintf(fout,"%d\n", NIND);
	fprintf(fout,"%d\n", numSNP);
	if (accmean(&AVE_numNIND) > 0)	fprintf(fout,"%f\n", 1.0/accmean(&AVE_numNIND));
	else		fprintf(fout,"-99\n");
	if ((numSNP > 0)&&(accmean(&pq) > 0))	fprintf(fout,"%f\n", f);
	else		fprintf(fout,"-99\n");
//	if (numSNP > 0)	fprintf(fout,"%f\n", fp);
//	else		fprintf(fout,"-99\n");

		for (b=1; b<=NBIN; b++)
		{
			if (n_[b] == 0)
			{
				if (b <= 5)	fprintf(fout,"0 0 0 %d 0 0\n", 2*b);
				else		fprintf(fout,"0 0 0 %d 0 0\n", ((NGEN/NBIN)*(b-5))+10);
			}
			else
			{
				if (b <= 5) fprintf(fout,"%lld %f %f %d %f %f\n", n_[b], cc_[b], d2_[b], 2*b, D2_[b], VV_[b]);
				else fprintf(fout,"%lld %f %f %d %f %f\n", n_[b], cc_[b], d2_[b], ((NGEN/NBIN)*(b-5))+10, D2_[b], VV_[b]);
			}
		}

	return(0);
}

/* ********************************************************************************************* */
