
/*GENLIB library routines */

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define FREEALL free_ivector(l3,1,m);free_ivector(l2,1,m);\
    free_ivector(l1,1,n+1);
#define NR_END 1
#define FREE_ARG char*
#define pi 3.14159265358979
#define infinity 999999
#define true 1
#define false 0
#define maxranges 5                    /* for numerical integration */
#define maxsimpsonpoints 1025          /* change libhdr if these changed */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define MBIG 1000000000
#define MSEED 161803398
#define FAC (1.0/MBIG)


typedef double (*fun1)();

FILE *seedfileptr, *tmoutfile, *fopen();


int inputcounter;
long seed;
long rseed;
int tracelevel;

long iy=0;
long iv[NTAB];

struct acc
{
   int n;
   double sum;
   double sumsq;
   int n2;
   double sum2;
   double sumsq2;
};

struct covacc
{
   int n;
   double sumxy, sumx, sumy, sumx2, sumy2;
};


/* fast uniform 0-1 generatot */
float ranqd2()
{
   unsigned long itemp;
   static unsigned long jflone = 0x3f800000;
   static unsigned long jflmsk = 0x007fffff;
   rseed = rseed*1664525L + 1013904223L;
   itemp = jflone | (jflmsk & rseed);
   return((*(float *)&itemp)-1.0);
}
   
double uniform()
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;
	
	if(seed<0||iff==0)
	{
		iff=1;
		mj=MSEED-(seed<0?-seed:seed);
		mj%=MBIG;
		ma[55]=mj;
		mk=1;
		for(i=1;i<=54;i++)
		{
			ii=(21*i)%55;
			ma[ii]=mk;
			mk=mj-mk;
			if(mk<0) mk+=MBIG;
			mj=ma[ii];
		}
		for(k=1;k<=4;k++)
			for(i=1;i<=55;i++)
			{
				ma[i] -=ma[1+(i+30)%55];
				if(ma[i]<0) ma[i]+=MBIG;
			}
		inext=0;
		inextp=31;
		seed=1;
	}
	if(++inext==56) inext=1;
	if(++inextp==56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if(mj<0) mj+=MBIG;
	ma[inext]=mj;
	seed=mj;
	return mj*FAC;
}






double ran2()
/* Returns a random number sampled from the uniform distribution
in range 0 to 1.  This is the ran2 function of p282 in
Numerical Recipes in C. */
{
	int j;
	long k;
	double temp;
	
	if(seed<=0||!iy)
	{
		if(-seed<1) seed=1;
		else seed= -seed;
		for(j=NTAB+7;j>=0;j--)
		{
			k=seed/IQ;
			seed=IA*(seed-k*IQ)-IR*k;
			if(seed<0) seed+=IM;
			if(j<NTAB) iv[j]=seed;
		}
		iy=iv[0];
	}
	k=seed/IQ;
	seed=IA*(seed-k*IQ)-IR*k;
	if(seed<0) seed+=IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j]=seed;
	if((temp=AM*iy)>RNMX) return RNMX;
	else return temp;
}

/*
double uniform()
A quick, but machine-dependent and dirty routine.
{
    return((double)rand()/65536/32768);
}
*/

#define ib1 1
#define ib2 2
#define ib5 16
#define ib18 131072

int irbit1(unsigned long *iseed)
/* routine to return random bits */
{
  unsigned long newbit;
  newbit = (*iseed &ib18) >> 17
  ^ (*iseed &ib5) >> 4
  ^ (*iseed &ib2) >> 1
  ^ (*iseed &ib1);
  *iseed = (*iseed<<1) | newbit;
  return (int)newbit;
}


getseed	()
{
   printf("Enter seed (-99 to read from file) ");
   manual: scanf("%ld", &seed);
   if (seed == -99)
   {
      seedfileptr = fopen("seedfile", "r");
      if (seedfileptr==0)
      {
         printf("No seedfile, enter seed please ");
         goto manual;
      }
      fscanf(seedfileptr, "%ld", &seed);
   }
   if (seed>0) seed = -seed;
   printf("Seed %ld\n", seed);
   rseed = seed; if (rseed<0) rseed = -rseed;   /* seed for ranqd2() */
/*   srand(seed);*/
}

getseedquick()
{
   seedfileptr = fopen("seedfile", "r");
   if (seedfileptr==0)
   {
      printf("No seedfile, enter seed please ");
      scanf("%ld", &seed);
   }
   else fscanf(seedfileptr, "%ld", &seed);
/*   srand(seed);*/
}


writeseed()
{
   double r, temp;
   seedfileptr = fopen("seedfile", "w");
   r = uniform();
   temp = floor(r*100000);
   seed = (long)temp;
   fprintf(seedfileptr, "%ld\n", seed);
}




initacc(a)
struct acc *a;
{
   a -> n = 0;
   a -> sum = 0;
   a -> sumsq = 0;
}


accum(struct acc *a, double x)
{
   a->n = a->n + 1;
   a->sum = a->sum + x;
   a->sumsq = a->sumsq + x*x;
}

double accmean(struct acc *a)
{
   return(a->sum/(double)(a->n));
}

double accsum(struct acc *a)
{
   return(a->sum);
}

double variance(struct acc *a)
{
   double num, denom;
   if (a->n == 0) return((double)-infinity);
   num = a->sumsq - (a->sum)*(a->sum)/((double)(a->n));
   denom = (double)(a->n) - (double)1;
   return(num/denom);
}


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! normal
! ------
! Returns random long real from the normal distribution.
!
! Parameters: mu = mean
!           sdev = standard distribution.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/



double normal(mu, sdev)
double mu, sdev;
{
  double u1, u2, r;
  u1= uniform();
  if (u1==0) u1 = 0.00001;            /*prevent fatal error*/
  if (u1==1) u1 = .999999;
  u2= uniform();
  r = sqrt (-2*log(u1)) * cos(2*pi*u2);
  return(r*sdev + mu);
}




/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! quick normal
! ------------
! Returns two uniform long reals from the normal distribution.
!
! Parameters: x1, x2 -> return values
!           mu = mean
!           sdev = standard distribution.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
quicknormal (x1, x2)
double *x1, *x2;
{
  double u1, u2, r;
  u1= uniform();
  if (u1==0) u1 = 0.00001;            /*prevent fatal error*/
  if (u1==1) u1 = .999999;
  u2 = uniform();
  *x1 = sqrt (-2.0*log(u1)) * cos(2.0*pi*u2);
  *x2 = sqrt (-2.0*log(u1)) * sin(2.0*pi*u2);
}



/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! wishart
! -------
!
! Returns two uniform real numbers from the wishart (bivariate) distribution.
! These are obtained by squaring correlated uniform numbers from the bivariate
! normal distribution.
!
! Parameters : z1 = return value for symmetrical distribution (i.e. metric)
!              z2 =  ,,     ,,    ,, negative sided distribution (i.e. fitness)
!              epsilon1 = sqrt(E(z1*z1))
!              epsilon2 = sqrt(E(z2*z2))
!              rho = correlation of absolute values.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
wishart(z1, z2, epsilon1, epsilon2, rho)
double *z1, *z2, epsilon1, epsilon2, rho;
{
   double sq;
   rho = sqrt(rho);
   sq = 1/sqrt(3.0);
   quicknormal(z1, z2);
   *z2 = *z1*rho + *z2*sqrt(1.0 - rho*rho);
   *z1 = sq*epsilon1*(*z1)*(*z1);
   *z2 = -sq*epsilon2*(*z2)*(*z2);
   if (uniform() > 0.5) *z1 = -(*z1);
}


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! bivnormal
! -------
!
! Returns two uniform real numbers from the bivariate normal distribution.
!
! Parameters : z1 = return value for symmetrical distribution (i.e. metric)
!              z2 = return value for symmetrical distribution (i.e. metric)
!              epsilon1 = sqrt(E(z1*z1))
!              epsilon2 = sqrt(E(z2*z2))
!              rho = correlation of absolute values.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
bivnormal(z1, z2, epsilon1, epsilon2, rho)
double *z1, *z2, epsilon1, epsilon2, rho;
{
   rho = sqrt(rho);
   quicknormal(z1, z2);
   *z2 = *z1*rho + *z2*sqrt(1.0 - rho*rho);
   *z1 = epsilon1*(*z1);
   *z2 = epsilon2*(*z2);
}



cubenormal(z1, z2, epsilon1, epsilon2, rho)
double *z1, *z2, epsilon1, epsilon2, rho;
{
   double sq;
   sq = 1.0/sqrt(16.0);
   rho = sqrt(rho);
   quicknormal(z1, z2);
   *z2 = *z1*rho + *z2*sqrt(1.0 - rho*rho);
   *z1 = sq*epsilon1*(*z1)*(*z1)*(*z1);
   *z2 = -sq*epsilon2*(*z2)*(*z2)*(*z2);
   if (*z2 > 0.0) *z2 = -*z2;
}



/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! doublegamma
! -----------
!
! Returns uniform long real from the double gamma distribution.
!
! Parameter: epsilon = sqrt(E(a*a))
!            proportion positive = proportion of distribution +ve
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
double doublegamma (epsilon, proppositive)
double epsilon, proppositive;
{
  double r, sigma;
  sigma = sqrt(epsilon*1.0/sqrt(3.0));
  r= normal(0.0, sigma);
  r= r*r;
  if (uniform() > proppositive) r= -r;       /*uniformly allocate sign*/
  return(r);
}


double gamdev(int ia, double epsilon)
/* returns random number from gamma distribution with integer parameter
taken from Numerical Recipes in C */
{
   int j;
   double x;
   if ((ia<1)||(ia>=6)) gabort("Parameter to gamdev() invalid ", (double)ia);
   x = 1.0;
   for (j=1; j<=ia; j++)
   {
      x *= uniform();
   }
   x = -log(x);
   x = (epsilon*x)/sqrt((double)(ia*(ia+1)));
   return x;
}









/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! init covaccumulation
! -----------------
!
! Sets to zero all variables in a mean and covariance accumulation record.
!
! Parameters:
!              a -> covariance accumulation record.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


initcovacc(a)
struct covacc *a;
{
   a->n = 0;
   a->sumx = 0;
   a->sumy = 0;
   a->sumxy = 0;
   a->sumx2 = 0;
   a->sumy2 = 0;
}


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! covaccumulate
! -------------
!
! Accumulates sums and sums of products.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

covaccum(a, x, y)
struct covacc *a;
double x, y;
{
   a->n = a->n + 1;
   a->sumx = a->sumx + x;
   a->sumy = a->sumy + y;
   a->sumxy = a->sumxy + x*y;
   a->sumx2 = a->sumx2 + x*x;
   a->sumy2 = a->sumy2 + y*y;
}



/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! covariance
! --------------------
!
!
! Returns the covariance calculated from an cov accumulation type record.
!
! Parameters: a-> accumulation record
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

double covariance(a)
struct covacc *a;
{
   if (a->n < 2)
   {
      return(-infinity);
   }
   return((a->sumxy - (a->sumx*a->sumy/(double)a->n))/((double)a->n - 1.0));
}



/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! correlation
! -----------
!
!
! Returns the correlation calculated from an cov accumulation type record.
!
! Parameters: a-> accumulation record
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

double correlation(a)
struct covacc *a;
{
   double num, denom, cov, v1, v2;
   if (a->n < 2)
   {
      return(-infinity);
   }
   cov = covariance(a);
   num = a->sumx2 - (a->sumx)*(a->sumx)/((double)(a->n));
   denom = (double)(a->n) - 1.0;
   v1 = num/denom;
   num = a->sumy2 - (a->sumy)*(a->sumy)/((double)(a->n));
   denom = (double)(a->n) - 1.0;
   v2 = num/denom;
   return(cov/sqrt(v1*v2));
}




/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! se
! --
!
! Computes the standard error in a record of type accum.
!
! Parameter: a -> record.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
double se(a)
struct acc *a;
{
   if (a->n < 1) return(-infinity);
   if (variance(a)/a->n < 0) return(-infinity);
   return(sqrt(variance(a)/a->n));
}



/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! print mse
! ---------
!
! Prints the mean and standard error in the record accumulation.
!
! Parameters: s -> string header.
!             r -> record>
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 
printmse(s, r)
char *s;
struct acc *r;
{

   printf("%s %f +/- %f\n", s, accmean(r), se(r));
}


gabort(s, r)
char *s;
int r;
{
   printf("ERROR %s %d\n", s, r);
   exit(1);
}



/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! print msd
! ---------
!
! Prints the mean and standard deviation.
!
! Parameters: s -> string header.
!             r -> record>
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 
printmsd(s, r)
char *s;
struct acc *r;
{

   printf("%s %f +/- %f\n", s, accmean(r), sqrt(variance(r)));
}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! generate poisson table
! ----------------------
!
! generates a table of probabilities for each whole number given a poisson
! distribution.
!
! Parametes: mean = mean of poisson.
!            last number = max. elements in table
!            table -> table of probabilities.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

generatepoissontable (mean, lastnumber, table, max)
double mean, *table;
int *lastnumber, max;
{
   int j;
   double prev;
   prev = exp(-mean);
   if (prev == 0) gabort("Insufficient resolution in generate poisson table", 0);
   table[0] = prev;
   for (j = 1; j <= max; j++)
   {
      *lastnumber = j;
      prev = prev*mean/(double)j;
      table[j] = table[j-1] + prev;
      if (1 - table[j] < 0.000001) break;
   }
}




/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! poisson
! -------
!
! Looks up a table of probabilities of the poisson distribution already
! set up.  Returns the index of this probability.  This discrete vaue
! will come from the posson distribution.
!
! Parameters: last number = last entry in the table
!             table      -> table of probabilities
!
! Returns   : integer from poisson distribution.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

int poisson (lastnumber, table)
int lastnumber;
double *table;
{
   int j;
   double u;
   u = uniform();
   for (j = 0; j<=lastnumber; j++)
   {
      if (u < table[j]) return(j);
   }
   return(lastnumber);
}


/*-------------------------------------------------------------------*/
/* binarysearchcumulative                                            */
/* ----------------------                                            */

/* Assumes array contains a cumulative density function of a discrete*/
/* distribution.  Returns an                                         */
/* integer randomly drawn from the density function.                 */

/* Parameters array -> cdf.                                          */
/*            size = no. of elements in array [1..size]              */
/*-------------------------------------------------------------------*/

int binarysearchcumulative(array, size)
double *array;
int size;
{
   int cur, pre, post;
   double r;
   r = uniform();
/*   printf("Uniform %f\n", r);*/
   pre = 1;
   post = size;
   cur = (size + 1)/2;
   do
   {
      if (array[cur] < r) pre = cur; else post = cur;
      cur = (pre + post)/2;
      if ((cur==post) || (cur==pre))
      {
         if (r < array[pre])
         {
/*            printf("Returning pre %d\n", pre);*/
            return(pre);
         }
         else
         {
/*            printf("Returning post %d\n", post);*/
            return(post);
         }
      }
   } while (size > 0);
}



/*----------------------------------------------------------------------*/
/*									*/
/* discrete 								*/
/* --------								*/
/*									*/
/* Returns random integer in range 1..n with equal probability.          */
/* Parameter: n = max. integer to return.				*/
/*									*/
/*----------------------------------------------------------------------*/
int discrete(n)
int n;
{
   return((int)(uniform()*(double)n) + 1);
}



   
/*----------------------------------------------------------------------*/
/*                            						*/
/* samplewithoutreplacement						*/
/* ------------------------						*/
/*									*/
/* Returns a random integer from those present in array, which is valid */
/* from 1 to limit.  Overwrites sampled integer with array[limit],      */
/* then decrements limit.						*/
/* Parameters: limit -> no. valid integers in array                     */
/*             array -> array[1..limit] of integers to sample           */
/* 									*/
/*----------------------------------------------------------------------*/
int samplewithoutreplacement(limit, array)
int *limit, *array;
{
   int index, res;
   index = discrete(*limit);
   res = array[index];
   array[index] = array[*limit];
   *limit = *limit - 1;
   return(res);
}
   
/* trap - asks for input to stop program */
trap()
{
   int i;
   printf("Enter any int to continue\n");
   scanf("%d", &i);
}


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! getint;
! -------
!
! Displays a prompt and reads in an integer.  Checks the range of the integer
! and stops the program if it is out of range.
!
! Parameters: s -> string prompt
!             i -> integer variable to input.
!           min =  max value
!           max = min value.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
getint(s, i, min, max)
char *s;
int *i, min, max;
{
   printf("%s", s);
   scanf("%d", i);
   if ((min != -infinity) && (*i < min))
   {
      printf("Value too small. Reading %s. Terminating.", s);
      printf("\n");
      printf("Min allowable: ");
      printf("%d", min);
      printf("\n");
      gabort("", 0);
   } 
   if ((max != infinity) && (*i > max)) 
   {
      printf("Value too large. Reading %s. Terminating.", s);
      printf("\n");
      printf("Max allowable: ");
      printf("%d", max);
      printf("\n");
      gabort("", 0);
   }
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! getint and skip;
! ----------------
!
! Displays a prompt and reads in an integer, skipping
! to the end of the line.  Checks the range of the integer
! and stops the program if it is out of range.
!
! Parameters: s -> string prompt
!             i -> integer variable to input.
!           min =  max value
!           max = min value.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
getintandskip(s, i, min, max)
char *s;
int *i, min, max;
{
   char ch;
   printf(s);
   printf(" ");
   scanf("%d", i);
   do
   {
     ch = getchar();
   }
   while (ch != '\n');
   if ((min != -infinity) && (*i < min))
   {
      printf("Value too small. Reading %s. Terminating.", s);
      printf("\n");
      printf("Min allowable: ");
      printf("%d", min);
      printf("\n");
      gabort("", 0);
   } 
   if ((max != infinity) && (*i > max)) 
   {
      printf("Value too large. Reading %s. Terminating.", s);
      printf("\n");
      printf("Max allowable: ");
      printf("%d", max);
      printf("\n");
      gabort("", 0);
   }
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
printreal(s,r,  x, y)
char *s;
double r;
int x, y;
{
   printf("%s %f\n", s, r);

}




/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! printint;
! ---------
!
! Prints out the integer with leading spaces x after the string s.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
printint (s, i, x)
char *s;
int i, x;
{
   printf("%s %d", s, i);
   printf("\n");
}


/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! acc2umulate
! ----------
!
! Accumulates the number of entries value, the sum entry and the sum of
! squares entry in an accumulation type record.
! Accumulated in the second accumulation counters.
!
! Parameters: a -> accumulation record.
!             x =  value to accumulate.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
acc2umulate (a, x)
struct acc *a;
double x;
{
   a->n2 = a->n2 + 1;
   a->sum2 = a->sum2 + x;
   a->sumsq2 = a->sumsq2 + x*x;
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! acc2umulated variance
! --------------------
!
!
! Returns the variance calculated from an accumulation type array.
! Uses second set of accumulation counters.
!
! Parameters: a-> accumulation record
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

double acc2umulatedvariance(a)
struct acc *a;
{
   if (a->n2 < 2 )
   {
      printf("In routine ACUUMULATE2\n");
      return(-infinity);
   }
   return((a->sumsq2 - (a->sum2*a->sum2/(double)a->n2))/((double)a->n2 - 1.0));
}




/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! se2
! --
!
! Computes the standard error in a record of type accum.
! Uses the second set of accumulation records.
!
! Parameter: a -> record.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
double se2(a)
struct acc *a;
{
   if (a->n2 < 1)
   {
      return(-infinity);
   }
   return(sqrt(acc2umulatedvariance(a)/(double)(a->n2)));
}



/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! acc2umulated mean
! ----------------
!
! Returns a mean value from an accumulation type record.
! Uses the second set of accumulation counters.
!
! Parameters: a ->> accumulation type record.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
double acc2umulatedmean(a)
struct acc *a;
{
   if (a->n2 == 0)
   {
      return(-infinity);
   }
   return(a->sum2/(double)a->n2);
}





/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! print mse2
! ---------
!
! Prints the mean and standard error in the record accumulation.
! If the mode is for multiple independent runs the accumulates into
! the second sest of accumulation counters.  Otherwise, or if the last
! run has been done, prints out the mean and s.e.
!
! Parameters: s -> string header.
!             r -> record>
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

printmse2(s, r, run, runs, multiplerunsflag)
char *s;
struct acc *r;
int run, runs, multiplerunsflag;
{
   if (multiplerunsflag == false )
   {
      printmse(s, r);
   }
   else
   {
      if (run == 1)
      {
          r->n2 = 0;
          r->sum2 = 0.0;
          r->sumsq2 = 0.0;
      }
      acc2umulate(r, accmean(r));
      r->sumsq = 0;
      r->sum = 0.0;
      r->n = 0.0;
      if (run == runs )
      {
         printf("%s %f", s, acc2umulatedmean(r));
         printf(" +/- ");
         printf("%f", se2(r));
         printf("*");
         printf("\n");
      }
   }
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! print tail
! ----------
!
! Prints a line of dashes.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

printtail()
{
   int j;
   printf("\n");
   for (j=1; j<=80; j++)
   {
      printf("-");
   }
   printf("\n\n");
}






/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! cointoss
! --------
!
! Uses random number to return false with probability 0.5 otherwise
! returns a non zero integer.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
int cointoss()
{
   if (uniform() > 0.5) return(0);
   return(1);
}



/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! printline
! ----------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
printline (s)
char *s;
{
   printf(s);
   printf("\n");
}

spaces(n)
int n;
{
   int i;
   for (i=1; i<=n; i++) printf(" ");
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! getreal and skip;
! ----------------
!
! Displays a prompt and reads in an real, skipping
! to the end of the line.  Checks the range of the integer
! and stops the program if it is out of range.
!
! Parameters: s -> string prompt
!             i -> integer variable to input.
!           min =  max value
!           max = min value.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
getrealandskip(s, i, min, max)
char *s;
double *i, min, max;
{
   char ch;
   printf(s);
   printf(" ");
   scanf("%lf", i);
   do
   {
     ch = getchar();
   }
   while (ch != '\n');
   if ((min != -infinity) && (*i < min))
   {
      printf("Value too small. Reading %s. Terminating.", s);
      printf("\n");
      printf("Min allowable: ");
      printf("%f", min);
      printf("\n");
      gabort("", 0);
   } 
   if ((max != infinity) && (*i > max)) 
   {
      printf("Value too large. Reading %s. Terminating.", s);
      printf("\n");
      printf("Max allowable: ");
      printf("%f", max);
      printf("\n");
      gabort("", 0);
   }
}

double calculaterealmean(a, n)
double *a;
int n;
{
   double sum;
   int i;
   sum = 0.0;
   for (i=1; i<=n; i++)
   {
      sum = sum + a[i];
   }
   return(sum/(double)n);
}
   
tracestart()
{
printf("Enter trace level ");
scanf("%d", &tracelevel);
}

trace(s, i)
{
   if (tracelevel==0) return(0);
   printf("%s %d\n", s, i);
}


outputrep(j, howoften)
int j, howoften;
{
  if ((double)j/(double)howoften  - floor((double)j/(double)howoften) == 0.0)
     printf("Completed iteration %d\n", j);
}


dummystart()
{
   inputcounter = 1;
}



dummyinput()
{
   inputcounter--;
   if (inputcounter <= 0)
   {
      printf("Enter int to continue ");
      scanf("%d", &inputcounter);
   }
}


double normalheight(x, mu, var)
double x, mu, var;
{
   double res;
   res  = (1.0/sqrt(2.0*pi*var))*exp(-(x-mu)*(x-mu)/(2*var));
   return(res);
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! simpson
! -------
!
! Computes the area under the function using simpsons rule for the set of
! points given.
!
! Parameters: x -> array of pairs of points.
!              points = no. of points.
!                   a = lower value
!                   v = upper value
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

double simpson (x, points, a, b)
double *x, a, b;
int points;
{
   int j, coeff;
   double  res, h;
   if (2*(points/2) == points)
   {
      printf("FATAL ERROR: Even no. of points for SIMPSON.\n");
      dummystart();
      dummyinput();
   }
   h = (b-a) / ((double)points-1.0);
   res = 0.0;
   for (j = 1; j<=points; j++)
   {
      if ((j == points) || (j==1 ))
      {
         coeff = 1;
      }
      else if (j == 2)
      {
         coeff = 4;
      }
      else if (coeff == 2) coeff = 4; else coeff = 2;
/*      %if trace level > 0 %then printstring("Coeff, x(j) ") %and write(coeff, 3) %and print(x(j), 2, 4) %and newline*/
      res = res + (double)coeff*x[j];
   }
   return((h/3.0) * res);
}




odd(i)
int i;
{
   int x;
   x = i/2;
   if ((double)(i)/2.0 - (double)x == 0.0) return(false); else return(true);
}




/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! solve quadratic
! ---------------
!
! Returns the roots of quadratic of parameters given.
!
! Parameters: a, b, c = parameters of quadratic.
!             root1, root2 = roots of quadratic.
!
! Returns : true => real roots.
!          false => no real root.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
int solvequadratic(a, b, c, r1, r2)
double a, b, c, *r1, *r2;
{
   double x;
   x = b*b - 4.0*a*c;
   if (x < 0.0) return(false);
   *r1 = (-b + sqrt(x))/(2.0*a);
   *r2 = (-b - sqrt(x))/(2.0*a);
   return(true);
}

   



skiptoendofline(fptr)
FILE *fptr;
{
   char ch;
   do
   {
      fscanf(fptr, "%c", &ch);
   }
   while (ch!='\n');
}



/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
! aitken
! ------
!
! Uses Aitken's formula to predict the asymptote from the three last points
! in the array given.
!
! Parameters: x-> array of points.
!               t = max. point in array.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
double aitken (a, t)
double *a;
int t;
{
   double x;
   x = a[t-2] - 2.0*a[t-1] + a[t];
   if (x == 0) return(a[t]);
   return(-((a[t-1] - a[t])*(a[t-1] - a[t]))/x + a[t]);
}

int factorial(i)
int i;
{
   int res;
   res = 1;   
   do
   {
      if (i==0) break;
      res = res*i;
      i--;
   }
   while (true==true);
   return(res);
}

double logfactorial(i)
int i;
{
   double res = 0.0;
   do
   {
      if (i==0) break;
      res = res + log((double)i);
      i--;
   }
   while (true==true);
   return(res);
}

double doublefactorial(x)
double x;
{
   double res;
   x = floor(x);
   res = 1.0;   
   do
   {
      if (x==0.0) break;
      res = res*x;
      x = x - 1.0;
   }
   while (true==true);
   return(res);
}


/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! gammaalpha
! ----------
!
! Returns the value of the parameter alpha of the gamma function.
!
! 
! Parameters: beta = parameter of gamma func.

!            epsilon = sqrt(E(a*a))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
double gammaalpha( beta, epsilon)
double beta, epsilon;
{
   double res;
   res = sqrt(beta*(beta + 1.0))/epsilon;
   return(res);
}




double asymptote_gamma(double z)
{
   double res;
   res = (z - 0.5)*log(z) -z +0.5*log(2.0*pi) + 1.0/(12.0*z) 
        - 1.0/(360.0*z*z*z) + 1.0/(1260.0*z*z*z*z*z) - 1.0/(1680.0*z*z*z*z*z*z*z);
   return(exp(res));
}
   
    
   
double series_gamma(double z)
/* return gamma(z) using series expansion 6.1.34 of Abramowitz and Stegun */
{
   #define maxc 26
   double c[maxc+1], res;
   int i;
   if (z > 2.0) return(asymptote_gamma(z));
   c[1] = 1.0;
   c[2] = 0.5772156649015329;
   c[3] =-0.6558780715202538;
   c[4] =-0.0420026350340952;
   c[5] = 0.1665386113822915;
   c[6] =-0.0421977345555443;
   c[7] =-0.0096219715278770;
   c[8] = 0.0072189432466630;
   c[9] =-0.0011651675918591;
   c[10]=-0.0002152416741149;
   c[11]= 0.0001280502823882;
   c[12]=-0.0000201348547807;
   c[13]=-0.0000012504934821;
   c[14]= 0.0000011330272320;
   c[15]=-0.0000002056338417;
   c[16]= 0.0000000061160950;
   c[17]= 0.0000000050020075;
   c[18]=-0.0000000011812746;
   c[19]= 0.0000000001043427;
   c[20]= 0.0000000000077823;
   c[21]=-0.0000000000036968;
   c[22]= 0.0000000000005100;
   c[23]=-0.0000000000000206;
   c[24]=-0.0000000000000054;
   c[25]= 0.0000000000000014;
   c[26]= 0.0000000000000001;
   res = 0.0;
   for (i=1; i<=26; i++)
   {
      res += c[i]*pow(z, (double)i);
/*      printf("gamma(z) %15.12f\n", 1.0/res); */
   }
   return(1.0/res);
}





/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! gammaht
! -------
!
! Returns the height at point a of a gamma distribution.
!
! Parameters: epsilon, beta = gamma parameters.
!             a = value to evaluate.
!
! Returns gamma fn.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
double gammaht(epsilon, beta, a)
double epsilon, beta, a;
{
   double gammabeta, res, a1, a2, a3, alpha;
   if (a==0.0)
   {
      printf("Illegal value (0.0) for gammaht()\n");
      dummyinput();
   }
   if (beta < 0.0) gabort("gammaht: invalid beta\n", beta);
   alpha = gammaalpha(beta, epsilon);
   gammabeta = series_gamma(beta);
   a1 = pow(alpha, beta);
   a2 = exp(-alpha*a);
   a3 = pow(a, beta - 1.0);
   res = a1*a2*a3/gammabeta;
   return(res);
}


/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! half points
! -----------
!
! This routine is called to check that the simpson numerical integration
! is converging satisfactorily.  It halves the number of points in
! the array specified.
!

! Parameters: a -> array of points.
!              points -> no. of points.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
halfpoints(a, points)
double *a;
int *points;
{
   int ptr1, ptr2, newpoints;
   if (odd(*points)==false)
   {
      printf("Half points: ERROR: no. of points for simpson must be odd.35\n");
      dummyinput();
   }
   ptr1 = 1;
   ptr2 = 1;
   do
   {
      a[ptr1] = a[ptr2];
      newpoints = ptr1;
      ptr1++;
      ptr2 = ptr2 + 2;
   }
   while (ptr2 <= *points);
   *points = newpoints;
}



setuppoints(pts, points, lower, upper, fn)
int pts;
double  *points, lower, upper;
fun1 fn;
{
      int i;
      double interval, x, a;
      interval = (upper - lower)/((double)pts - 1.0);
      x = lower;
      for (i = 1; i<=pts; i++)
      {
         a = fn(x);
         points[i] = a;
         x = x + interval;
      }

}


/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
! numerical integ
! ---------------
!
! Numerically integrates the function fn and returns the result in res.
! Integrates over a number of ranges with different numbers of points per
! range.  Points must be a number like 65, 129, 257, 513... in each
! range.  Halves number of points to do integration 3 times as a check on
! convergence
!
! Parameters:
!
!
! start, finish: values limiting the function
! fn(x): is the function to integrate
! ranges: no. of ranges of points to evaluate.
! points per range: no of points in each range.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

double numericalinteg(start, finish, fn, ranges, pointsperrange)
  double start, finish;
  fun1 fn;
  int ranges, *pointsperrange;
{
   double lower, upper, r;
   double res[maxranges+2][4];
/*        The first dimension is contains the integral for each range, the
       last element containing the total.  The second dimension contains
       results for each halving of the points.
*/
   double points[maxsimpsonpoints+1];
   int i, j, k, fac, p;
   static int tmoutfileflag = 0;           
   if (tmoutfileflag==0)
   {
      tmoutfile = fopen("numericalinteg.out", "w");
      tmoutfileflag = 1;
   }
/*   printf("Start %f, finish %f, ranges %d\n", start, finish, ranges);*/
   upper = start;
   for (i = 1;  i<= ranges; i++)
   {
      lower = upper;
      upper = lower + (finish - start)/(double)ranges;
      p = pointsperrange[i];
      if (p > maxsimpsonpoints) gabort("too many points for simpson",p);
      setuppoints(p, points, lower, upper, fn);
      for (j = 1; j<=3; j++)
      {
         r = simpson(points, p, lower, upper);
         res[i][j] = r;
         if (j != 3) halfpoints(points, &p);
      }
   }
   for (i = 1; i<=3; i++)
   {
      res[ranges+1][i] = 0.0;
      for (j = 1; j<=ranges; j++)
         res[ranges+1][i] = res[ranges+1][i] + res[j][i];
   }
   fprintf(tmoutfile, "Results of numerical integration\n");
   for (i=1; i<=ranges; i++) fprintf(tmoutfile, " Pts   Range %d", i);
   fprintf(tmoutfile, "     All");
   fprintf(tmoutfile, "\n");
   fac = 1;
   for (i = 1; i<=3; i++)
   {
      for (j = 1; j<=ranges+1; j++)
      {
         if (j!=ranges+1) fprintf(tmoutfile, "%4d", (pointsperrange[j]-1)/fac + 1);
         fprintf(tmoutfile, "%10.4f", res[j][i]);
      }
      fac = fac*2;
      fprintf(tmoutfile, "\n");
   }
   return(res[ranges+1][1]);
}



FILE *openforread(char *str)
{
   FILE *f;
   f = fopen(str, "r");
   if (f==0)
   {
      printf("ERROR: File %s not found.\n", str);
      exit(0);
   }
   else printf("Opened file %s for read.\n", str);
   return(f);
}



 
printmsefile(f, s, r)
FILE *f;
char *s;
struct acc *r;
{

   fprintf(f, "%s %f +/- %f\n", s, accmean(r), se(r));
}




quadratic_regression(double n, double sumx, double sumx2, double sumy, double sumx3,
  double sumxy, double sumx4, double sumx2y, double *b1, double *b2, double *b3)
{
/* Solve system of 3 simultaneous equations:

a*b1 + b*b2 + c*b3 = d         (1)
e*b1 + f*b2 + g*b3 = h         (2)
i*b1 + j*b2 + k*b3 = l         (3)

The actual parameters are appropriate for solving a quadratic regression -

y = b3*x^2 + b2*x + b3

See assignments below.

*/
   double a, b, c, d, e, f, g, h, i, j, k, l;
   double z1, z2, z3, z4, z5, z6;
   a = n;
   b = sumx;
   c = sumx2;
   d = sumy;
   e = b;
   f = c;
   g = sumx3;
   h = sumxy;
   i = c;
   j = g;
   k = sumx4;
   l = sumx2y;
/*   printf("\n%lf %lf %lf %lf\n", a, b, c, d);
   printf("%lf %lf %lf %lf\n", e, f, g, h);
   printf("%lf %lf %lf %lf\n\n", i, j, k, l);
*/

/* (1) - (2) -> b1*z1 + b2*z2 = z3  (4)     */
   z1 = a/c - e/g;
   z2 = b/c - f/g;
   z3 = d/c - h/g;

/* (1) - (3) -> b1*z4 + b2*z5 = z6  (5)     */
   z4 = a/c - i/k;
   z5 = b/c - j/k;
   z6 = d/c - l/k;


/* (4) - (5) */

   *b1 = (z3/z2 - z6/z5)/(z1/z2 - z4/z5);

   *b3 = ((d - a*(*b1))/b - (h - e**b1)/f) / (c/b - g/f);
  
   *b2 = (d - c*(*b3) - a*(*b1))/b;
}



int findtext(char *s, FILE *inptr)
{
   int len, i, curchar, dum;
   char c;
   curchar = 0;
   len = 0;
   for (i=0; i<=100; i++)
   {
      if (s[i] == '\0') break;
      len++;
   }
   for (;;)
   {
      c = getc(inptr);
      if (c == EOF) return(false);
      if (c==s[curchar])
      {
         curchar++;
         if (curchar==len) return(true);
      }
      else curchar = 0;
   }
   return(true);
}


/* ***************************************************************************** */
/* ***************************************************************************** */
/* *************** nrutil.c **************************************************** */
/* ***************************************************************************** */
/* ***************************************************************************** */


void nrerror(char error_txt[])
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_txt);
	fprintf(stderr,"...now exiting to system\n");
	exit(1);
}

float *vector(long nl, long nh)
{
	float *v;
	
	v=(float *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(float)));
	if(!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
{
	int *v;
	
	v=(int *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(int)));
	if(!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
{
	unsigned char *v;
	
	v=(unsigned char *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if(!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
{
	unsigned long *v;
	
	v=(unsigned long *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(unsigned long)));
	if(!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
{
	double *v;
	
	v=(double *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(double)));
	if(!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(long nrl,long nrh,long ncl,long nch)
{
	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	
	m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if(!m) nrerror("allocation failure 1 in matrix()");
	m+=NR_END;
	m-=nrl;
	m[nrl]=(float *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if(!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl]+=NR_END;
	m[nrl]-=ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	return(m);
}

double **dmatrix(long nrl,long nrh,long ncl,long nch)
{
	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	
	m=(double **)malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if(!m) nrerror("allocation failure 1 in dmatrix()");
	m+=NR_END;
	m-=nrl;
	m[nrl]=(double *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if(!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
	m[nrl]+=NR_END;
	m[nrl]-=ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	return(m);
}

int **imatrix(long nrl,long nrh,long ncl,long nch)
{
	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
	
	m=(int **)malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if(!m) nrerror("allocation failure 1 in imatrix()");
	m+=NR_END;
	m-=nrl;
	m[nrl]=(int *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if(!m[nrl]) nrerror("allocation failure 2 in imatrix()");
	m[nrl]+=NR_END;
	m[nrl]-=ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	return(m);
}

float **submatrix(float **a,long oldrl,long oldrh,long oldcl,long oldch,long newrl,long newcl)
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldch-oldcl+1;
	float **m;
	
	m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if(!m) nrerror("allocation failure in submatrix()");
	m+=NR_END;
	m-=newrl;
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
	return(m);
}

float **convert_matrix(float *a,long nrl,long nrh,long ncl,long nch)
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	
	m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if(!m) nrerror("allocation failure in convert_matrix()");
	m+=NR_END;
	m-=nrl;
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<=nrow;i++,j++) m[j]=m[i-1]+ncol;
	return(m);
}

float ***f3tensor(long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;
	
	t=(float ***)malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if(!t) nrerror("allocation failure 1 in f3tensor()");
	t+=NR_END;
	t-=nrl;
	t[nrl]=(float **)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if(!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl]+=NR_END;
	t[nrl]-=ncl;
	t[nrl][ncl]=(float *)malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if(!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl]+=NR_END;
	t[nrl][ncl]-=ndl;
	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++)
	{
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol+ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	return(t);
}

void free_vector(float *v,long nl,long nh)
{
	free((FREE_ARG)(v+nl-NR_END));
}

void free_ivector(int *v,long nl,long nh)
{
	free((FREE_ARG)(v+nl-NR_END));
}

void free_cvector(unsigned char *v,long nl,long nh)
{
	free((FREE_ARG)(v+nl-NR_END));
}

void free_lvector(unsigned long *v,long nl,long nh)
{
	free((FREE_ARG)(v+nl-NR_END));
}

void free_dvector(double *v,long nl,long nh)
{
	free((FREE_ARG)(v+nl-NR_END));
}

void free_matrix(float **m,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG)(m[nrl]+ncl-NR_END));
	free((FREE_ARG)(m+nrl-NR_END));
}

void free_dmatrix(double **m,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG)(m[nrl]+ncl-NR_END));
	free((FREE_ARG)(m+nrl-NR_END));
}

void free_imatrix(int **m,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG)(m[nrl]+ncl-NR_END));
	free((FREE_ARG)(m+nrl-NR_END));
}

void free_submatrix(float **b,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG)(b+nrl-NR_END));
}

void free_convert_matrix(float **b,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG)(b+nrl-NR_END));
}

void free_f3tensor(float ***t,long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
	free((FREE_ARG)(t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG)(t[nrl]+ncl-NR_END));
	free((FREE_ARG)(t+nrl-NR_END));
}


/* ***************************************************************************** */
/* ***************************************************************************** */
/* ******* Linear Programming routine ****************************************** */
/* ***************************************************************************** */
/* ***************************************************************************** */


void simplx(float **a, int m, int n, int m1, int m2, int m3, int *icase,
    int izrov[], int iposv[])
{    
    void simp1(float **a, int mm, int ll[], int nll, int iabf, int *kp,
        float *bmax);
    void simp2(float **a, int n, int l2[], int nl2, int *ip, int kp, float *q1);
    void simp3(float **a, int i1, int k1, int ip, int kp);
    int i,ip,ir,is,k,kh,kp,m12,nl1,nl2;
    int *l1,*l2,*l3;
    float q1,bmax;

    if (m!= (m1+m2+m3)) nrerror("Bad input constraint counts in simplx");
    l1=ivector(1,n+1);
    l2=ivector(1,m);
    l3=ivector(1,m);
    nl1=n;
    for (k=1; k<=n;k++) l1[k]=izrov[k]=k;
    nl2=m;
    for (i=1; i<=m;i++) {
        if (a[i+1][1] < 0.0) nrerror ("Bad input tableau in simplx");
        l2[i]=i;
        iposv[i]=n+i;
    }
    for (i=1; i<=m2;i++) l3[i]=1;
    ir=0;

    if (m2+m3) {
        ir=1;

        for (k=1; k<=(n+1); k++) {
	    q1=0.0;
	    for (i=m1+1;i<=m;i++) q1 += a[i+1][k];
	    a[m+2][k] = -q1;
        }
        do {
	    simp1(a,m+1,l1,nl1,0,&kp,&bmax);
	    if (bmax <= EPS && a[m+2][1] < -EPS) {
	        *icase = -1;

	        FREEALL return;
           } else if (bmax <= EPS && a[m+2][1] <= EPS) {
	       m12=m1+m2+1;

	       if (m12 <= m) {
	           for (ip=m12;ip<=m;ip++) {
		       if (iposv[ip] == (ip+n)) {
		           simp1(a,ip,l1,
			       nl1,1,&kp,&bmax);
		           if (bmax > 0.0)
			       goto one;
		       }
	           }
	       }
	       ir=0;
	       --m12;
	       if (m1+1 <= m12)
	           for (i=m1+1;i<=m12;i++)
		       if (l3[i-m1] == 1)
		          for (k=1;k<=n+1;k++)
	    		       a[i+1][k] = -a[i+1][k];
	       break;
           }
           simp2(a,n,l2,nl2,&ip,kp,&q1);
           if (ip == 0) {
	       *icase = -1;
	       FREEALL return;
           }
     one:  simp3(a,m+1,n,ip,kp);
           if (iposv[ip] >= (n+m1+m2+1)) {
	       for (k=1;k<=nl1;k++)
	           if (l1[k] == kp) break;
	       --nl1;
	       for (is=k;is<=nl1;is++) l1[is]=l1[is+1];
	       ++a[m+2][kp+1];
	       for (i=1;i<=m+2;i++) a[i][kp+1] = -a[i][kp+1];
           } else {
	       if (iposv[ip] >= (n+m1+1)) {
	           kh=iposv[ip]-m1-n;
	           if (l3[kh]) {
	 	       l3[kh]=0;
		       ++a[m+2][kp+1];
		       for (i=1;i<=m+2;i++)
		           a[i][kp+1] = -a[i][kp+1];
	           }
	       }
           }
           is=izrov[kp];
           izrov[kp]=iposv[ip];
           iposv[ip]=is;
        } while (ir);
    }

    for (;;) {
        simp1(a,0,l1,nl1,0,&kp,&bmax);
        if (bmax <= 0.0) {
	    *icase=0;
	    FREEALL return;
        }
        simp2(a,n,l2,nl2,&ip,kp,&q1);
        if (ip == 0) {
	    *icase=1;
	    FREEALL return;
        }
        simp3(a,m,n,ip,kp);
        is=izrov[kp];
        izrov[kp]=iposv[ip];
        iposv[ip]=is;
    }
}



void simp1(float **a, int mm, int ll[], int nll, int iabf, int *kp,
    float *bmax)
{
    int k;
    float test;

    *kp=ll[1];
    *bmax=a[mm+1][*kp+1];
    for (k=2;k<=nll;k++) {
	if (iabf == 0)
	    test=a[mm+1][ll[k]+1]-(*bmax);
	else
	    test=fabs(a[mm+1][ll[k]+1])-fabs(*bmax);
	if (test > 0.0) {
	    *bmax=a[mm+1][ll[k]+1];
	    *kp=ll[k];
	}
    }
}



void simp2(float **a, int n, int l2[], int nl2, int *ip, int kp, float *q1)
{
    int k,ii,i;
    float qp,q0,q;

    *ip=0;
    for (i=1;i<=nl2;i++) {
	if (a[l2[i]+1][kp+1] < -EPS) {
	    *q1 = -a[l2[i]+1][1]/a[l2[i]+1][kp+1];
	    *ip=l2[i];
	    for (i=i+1;i<=nl2;i++) {
		ii=l2[i];
		if (a[ii+1][kp+1] < -EPS) {
		    q = -a[ii+1][1]/a[ii+1][kp+1];
		    if (q < *q1) {
			*ip=ii;
			*q1=q;
		    } else if (q == *q1) {
			for (k=1;k<=n;k++) {
			    qp = -a[*ip+1][k+1]/a[*ip+1][kp+1];
			    q0 = -a[ii+1][k+1]/a[ii+1][kp+1];
			    if (q0 != qp) break;
			}
			if (q0 < qp) *ip=ii;
		    }
		}
	    }
	}
    }
}



void simp3(float **a, int i1, int k1, int ip, int kp)
{
    int kk,ii;
    float piv;

    piv=1.0/a[ip+1][kp+1];
    for (ii=1;ii<=i1+1;ii++)
	if (ii-1 != ip) {
	    a[ii][kp+1] *= piv;
	    for (kk=1;kk<=k1+1;kk++)
		if (kk-1 != kp)
		    a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1];
	}
    for (kk=1;kk<=k1+1;kk++)
	if (kk-1 != kp) a[ip+1][kk] *= -piv;
    a[ip+1][kp+1]=piv;
}



