/*
"fsimII.04.c"
FsimII
-------
Program created :	05/15/03	Anoop John
Last modified	:	07/21/03	Anoop John

Pseudocode
----------
{
	Read Control File
	Initialize data and counts
	Run Algorithm for initRunGen generations
		Repeat for repNum replicates
		{
			Run Algorithm for preRunGen generations
				Run Algorithm for repRunGen generations 
					Collect Data from replicate and output to file
					}
}
Algorithm
---------
{
	doMutations
	{
		For all individuals
		{
			For 1, 2, 3 codon positions
			{
				1->0 mutations	( ! remember to make same changes in both sections)
				0->1 mutations	( ! remember to make same changes in both sections)
			}
		}
	}
	doRecombination
	{
		For all pairs
		{
			find and sort positions for recombination
				swap alternate segments
			{
				if polymorphic site
					then if different
						then exchange data
						}
		}
	}
	getNextGenSeqs
	{
		Find Fitness of individuals
		Do multinomial sampling
		create next generation
	}
	updateCounts
	{
		Update counts of data being tracked
	}
}
*/
/*---------------------------------------------------------------------------------------*/
#define DEBUGMODE
#undef DEBUGMODE
/**/
#define VERBOSE
#undef VERBOSE
/**/
/*---------------------------------------------------------------------------------------*/
/*macros*/
/*options*/
#define DO_WRITE		1				
#define DONT_WRITE		0
#define DO_RESET		1	/* init data after prerun - 0 used for some pop fluct models */
#define DONT_RESET		0
#define INIT_RUN		0
#define PRE_RUN			1
#define REP_RUN			2
/*parameters*/
#define POP_MODEL		"constPopSize"		/*population fluctuation model */
#define MULTIHITFACTOR	10					/*every MULTIHITFACTOR*hitnum iterations check to see if there are free sites*/
#define INITMCU			0.6					/*Initial MCU value used to initialize the sequences                                  CHANGE THIS TO EQUILIBRIUM VALUE FOR MCU MODEL */
#define S_MODEL			"sampleAtEveryMutation"	/*sampling model for s values*/

/*constants*/
#define MAXSEQLEN		1500				/*maximum length of sequence has to be multiple of 3*/
#define MAXSAMPLENUM	25					/*maximum samplesize*/	
#define MAXSEQNUM		5000				/*maximum number of sequences even when fluctuating*/
#define MINSEQNUM		500					/*minimum number of sequences even under population fluctuation*/
#define MAXGEN			500000				/*maximum of inirungen, prerungen and reprungen*/
#define MAXFIXATIONS	200000				/*maximum number of fixations permitted in a single initrun or prerun or reprun*/
#define FLOATPRECISION	0.000001			/*precision to be used for calculations*/
#define PI				3.1415926535897932384626433832795
#define MAXLINECHAR		10000				/*maximum characters in a line*/
#define	MAXFILENAME		301					/*maximum length of filename*/
#define IDUM_INIT		-123341				/*value used to initialize the random number generator*/

/*---------------------------------------------------------------------------------------*/
/*structure definitions*/
struct seq_str
{
	char *site;					/*string of bases for each sequence*/
	long baseCount0[3+1];			/*baseCount0[j] is the number of 0s at sitePos where (sitePos+2)%3+1=j*/
	long baseCount1[3+1];			/*baseCount1[j] is the number of 1s at sitePos where (sitePos+2)%3+1=j*/
	long tot_1_ct;
};
struct siteData_str
{	
	char ancSite;							/*ancestral base at current site*/
	long freqs0;					/*frequency of 0s in the population at current site*/
	long freqs1;					/*frequency of 1s in the population at current site*/
	long fixNum01;					/*number of fixations from 0 to 1 at current site in the population during current replicate */
	long fixNum10;					/*number of fixations from 1 to 0 at current site in the population during current replicate */
	long mutNum01;					/*number of mutations from 0 to 1 at current site in the population during current replicate */
	long mutNum10;					/*number of mutations from 1 to 0 at current site in the population during current replicate */
	long recNum;					/*number of crossovers at current site during current replicate */
};
double *w_dat_static;
/*---------------------------------------------------------------------------------------*/
/*debug tools*/
#define pf(a) 		printf("%.20f\t", a)
#define pl(a) 		printf("%ld\t", a)
#define ps(a) 		printf("%s\t", a)
#define pa(a) 		printf("%c\t", a)
#define tab			printf("\t")
#define ent			printf("\n")
#define help		printf(" O'er h're")
#define halt		do{printf("\nPress enter to continue!\n");getchar();}while(0)
// #undef halt 
// #define halt		do{printf("\nPress any key to continue!\n");}while(0)

/*---------------------------------------------------------------------------------------*/
/*Error function implemented as macro*/
#define errorOut(a)	do{printf("\nError! ");printf a;printf("\nexiting program.....");getchar();exit(1);}while(0)
/*---------------------------------------------------------------------------------------*/
/*macro function to find codon position given site position*/
#define sitePosToCodPos(a)	((((a)+2)%3)+1)
/*---------------------------------------------------------------------------------------*/
/*global variables*/

/*flags read from the control file*/
long noMultiHit;								/*1=no multiple hits 0=permit multiple hits*/

/*parameters read from the control file*/
long repNum;						/*Total number of replicates for the simulation*/

long initRunGen;					/*number of repRunGen for the initial run */
long preRunGen;					/*number of repRunGen for the prerun before each reprun*/
long repRunGen;					/*number of generations for the replicate run*/	

long initSeqNum;					/*Initial number of sequences*/
long seqLen;						/*Length of sequences as number of bases*/
long sampleNum;					/*Number of sequences to be sampled to collect data*/

double c;									/*per site per generation recombination rate*/
double u10;									/*the 1->0 mutation rate common for all positions*/
double u01;									/*the 0->1 mutation rate common for all positions*/
double Nes2;								/* 2Nes for MCU model */

long idumVal, idum_init;						/*default ran2 seed. if seed given!=0 then it will be used*/	
long *idum=&idumVal;						
char outFileName[MAXLINECHAR+1], outFileName2[MAXLINECHAR+1], outFile_segsites[MAXLINECHAR+1], outFile_MCU[MAXLINECHAR+1];			/*name of output file*/
char outFileExtn[MAXLINECHAR+1];			/*outputfilename extension*/
char ctlFileName[MAXLINECHAR+1];			/*name of control file*/
char version[11];							/*version string will be appended to filename before extension*/
char fsimVersion[11];						/*version of fsim. the above is the version of the simulation*/
char comments[MAXLINECHAR+1];				/*comments - will be output into the outputfile*/
char machineName[MAXLINECHAR+1];			/*machine on which the program is running*/
FILE *fpOut, *fpOut2, *fpOut_segsites, *fpOut_MCU;								/*output file pointer*/

/*variables used for the simulation - data being processed*/
struct seq_str *seq, *newSeq;
struct siteData_str	siteData[MAXSEQLEN+1];	/*array of sitedata structures */
double sArray[MAXSEQLEN+1];					/*array of s values for the sequence*/

long preRunPopFluct[MAXGEN+1];		/*array to store the population sizes during each generation in the prerun*/
long popSizeArray[MAXGEN+1];		/*array used to store the population sizes for each gen for a given cycle*/
long curPopSize;					/*The current population size at any time*/
long curGenNum;					/*the current generation number*/
long curRepNum;					/*the current replicate*/
long sampleSeqs[MAXSAMPLENUM+1];	/*indices of sequences that are to be sampled for the current replicate*/		

/*arrays of flags for processing data*/
char curSeqHitSite[MAXSEQLEN+1]; 			/*1=site already hit in current mutagenesis of a gene, 0=not*/
char curGenHitSite[MAXSEQLEN+1];			/*1=site already hit in current generation (across genes) 0=not*/
char polySite[MAXSEQLEN+1];					/*1=site is segregating in population, 0=not*/		
char crossHitSite[MAXSEQLEN+1];				/*1=site already selected as crossover site in a given generation, 0=not*/

/*variables for data collection from each replicate. */
long sampleFix01[3+1];				/*Fixations of 0->1 due to sampling (fixed in sample, but not in pop) */
long sampleFix10[3+1];				/*Fixations of 1->0 due to sampling*/
long sampleFixNum;					/*Total number of fixations due to sampling */
long sampleFixPos[MAXSEQLEN+1];	/*Position of fixed (due to sampling) sites in the sample*/
long sampleSeg01[3+1];				/*segregating sites in the sample with ancSite 0*/
long sampleSeg10[3+1];				/*segregating sites in the sample with ancSite 1*/
long sampleSegNum;					/*Total number of segregating sites in the sample*/
long sampleSegPos[MAXSEQLEN+1];	/*Position of segregating sites in the sample*/

long popRecNum;					/*Total number of crossovers in the population across all sites */

/*for every fixation in the population, store the derived state, the position in the sequence, and the fitness effect */
long popFixNum;					/*number of mutations fixed in the replicate - used to access the following arrays*/
char fixSite[MAXFIXATIONS+1];				/*the site that got fixed ie. if 0 was fixed then 0 else 1*/
long fixSitePos[MAXFIXATIONS+1];	/*Positions of fixed sites in current replicate*/
double fixSiteSval[MAXFIXATIONS+1];			/*S values of sites fixed in current replicate*/

long repTime;						/*time in seconds for the current replicate*/

long genMultiHitCount;				/*Total of multiple hits in a generation across all generations for current rep*/
long segMultiHitCount;				/*Total of multiple hits at segregating sites across all generations*/
/*---------------------------------------------------------------------------------------*/
/*integer random number generator*/
#define getRandLong_(min, max)	((min)+(long)((1.0-ran2())*((max)-(min)+1)))	/* careful about which ran2 is being used !!!!!!!!!!!!! */
/*---------------------------------------------------------------------------------------*/
/*function declarations*/
long getCtlFile(void);
long checkParameters(void);
long writeFsimInfo(void);
long initData(void);
long doFsim(long genNum, long flagReset, long repType);
long getPopScenario(long repType, long genNum);
long resetData(void);
long doMutations(void);
long getHitSites(char origSite, long baseCount, long curSeq, long codPos, long* hitSites, long numHits);
long getPosNthChar(char *site, long codPos, long hitCod, char origSite);
long doRecombination(void);
int compareUL(const void*a, const void*b);
long getNextGenSeqs(void);

long updateCounts(void);
double sampleSval(long sitePos, char ancSite);
long writeRepData(void);
long getRepData(void);
long getSample(void);
long endSim(void);
long checkBaseCounts(long pos);

double ran2(void);
double poidev(double xm);
double gammln(double xx);
long getTime(long);
void multdev(double inn[], long k, long n, long nn[]);
double bnldev(double pp, long n);

FILE* fileOpen(char *path, char *mode);
char* getFileString(FILE *fp, char *string, long maxStringLength);
long readDataWithComments(FILE* fp, char*format, void* target);
void* memAlloc(size_t nmemb, size_t membsize, char* membname);

void get_static_w_dat(void);
long getNextGenSeqs_MCU(void);

/*---------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*safe memory allocation function. tries to allocate memory and if error then 
 print error with variable name and exit
 nmemb - the number of elements to be allocated
 memmsize - size of each element
 membname - name of the variable for which memory is being allocated for errorout
 */
void* memAlloc(size_t nmemb, size_t membsize, char* membname)
{
	void* ptr;
	static long totmem=0;
	if(!nmemb&&!membsize)
		return &totmem;
	ptr=calloc(nmemb, membsize);
	if(ptr)
	{
		totmem+=nmemb*membsize;
		return ptr;
	}
	else
		errorOut(("Unable to allocate %ld bytes for %s", nmemb*membsize, membname));
	return ptr;
}
/*--------------------------------------------------------------------------------*/
/*Safe file open function. can check for OS specific length issues with filenames
 for mac the length is 31 and for OSX it is 255
 path - the path to the file+filename or just filename if in same folder
 mode - the fopen modes
 */
FILE* fileOpen(char *path, char *mode)
{
	FILE *fp;
	long i, j;
	
	i=strlen(path);
	j=0;
	while(i>0)
	{
		i--;
		if(path[i]==':'||path[i]=="\\"[0]||path[i]=='/')
			break;
		j++;
	}
	if(j>MAXFILENAME)
		errorOut(("Filename %s exceeded the maximum limit set by the operating system", &path[i]));
	fp=fopen(path, mode);
	if(fp)
		return fp;
	else
		errorOut(("Unable to open %s in %s mode", path, mode));
	return fp;
}
/*--------------------------------------------------------------------------------*/
/*function to safely read a string with spaces from a file. checks for the maximum
 size of the string that can be read into the variable passed
 fp - file pointer at position where string is to be read
 string - the pointer to the location where the string is to be stored
 maxStringLength - the maximum length of the string that can be read into the variable
 created 		:06/06/03	- Anoop John 
 last modified	:06/06/03	- Anoop John
 */
char* getFileString(FILE *fp, char *string, long maxStringLength)
{
	long curLength;
	curLength=0;
	
	/*check for null pointer*/
	if(string==0)
		errorOut(("Null pointer passed to getFileString"));
	
	/*try to read till max length or end of file*/
	while(curLength<maxStringLength&&!feof(fp))
	{
		/*read char by character*/
		string[curLength]=fgetc(fp);
		/*read till end of line is read. break at end of line*/
		if(string[curLength]=='\n'||string[curLength]=='\r')
		{
			string[curLength]=0;
			break;
		}
		else
		{
			/*if new line is not reached after max length has been read then error*/
			curLength++;
			if(curLength==maxStringLength)
				errorOut(("String Length exceeded maximum allowed(%lu) in getFileString", maxStringLength));
		}
	}
	/*if end of file is reached error*/
	if(feof(fp))
		errorOut(("Unexpected end of file in getFileString"));
	return string;
}
/*--------------------------------------------------------------------------------*/
/*Function that reads in one data item from a file given the format
 as a string like in the normal scanf function
 fp 		- the pointer to the file from which it should read data from 
 format 	- equivalent to scanf format but only one variable will be read
 target  - address of the variable into which the data will be read into
 created 		:02/01/03	- Anoop John 
 last modified	:02/01/03	- Anoop John
 */
long readDataWithComments(FILE* fp, char*format, void* target)
{
	char tempc, foundeof=0;
	
	/*reads the data*/
	if(fscanf(fp, format, target)!=1)
		errorOut(("Corrupted input file or invalid format string(%s)", format));
	
	/*skips till end of line either CR or LF*/
	tempc=0;
	while((tempc!='\n')&&(tempc!='\r'))
	{
		tempc=fgetc(fp);
		if(foundeof)		/*Added check for ending CR or LF 01-27-03*/
			errorOut(("Missing ending CR or LF character in last line.\nPossible corrupted file "));
		if(feof(fp))
			foundeof=1;
	}
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Multinomial deviate function. Function fills in the array nn with the number of 
 successes for the corresponding probability of successes given in the inn array and
 the total number of experiments n. k is the size of the two arrays.(1 based arrays)
 Function Created:  Hiroshi Akashi
 */
void multdev(double inn[], long k, long n, long nn[])
{
	long i, r;
	double sum;
	sum = 0.0;
	for(i=1;i<=k;i++)
	{
		if( (n > 0) && (sum < 1.0) )
		{
			r = (long)(bnldev(inn[i]/(1.0 - sum), n));
			nn[i] = r;
			sum += inn[i];
			n -= r;
		}
		else
		{
			nn[i] = 0;
		}
	}
	return;
}
/*--------------------------------------------------------------------------------*/
/*Function returns the binomial deviate. Given the probability pp and the number of
 experiments the function returns the number of successes(long) as a double value.
 Function taken from Numerical Recipes for C IIed.
 */
double bnldev(double pp, long n)
{
	long j;
	static long nold=(-1);
	double am, em, g, angle, p, bnl, sq, t, y;
	static double pold=(-1.0), pc, plog, pclog, en, oldg;
	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (ran2() < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran2();
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran2();
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=(long)(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
										  -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran2() > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}
/*--------------------------------------------------------------------------------*/
/*Modified Ran2 function. Original function from the Numerical Recipes in C IIed. 
 This function uses a global idum. returns random number in the range (0, 1].
 */
#define IM1 2147483563
#define IM2 2147483399
#define AM2 (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 256
#define NDIV2 (1+IMM1/NTAB)
/*
 uncomment these lines as well as the ones at the end
 of the function to get (0, 1) output.
 */
/*
 #define EPS DBL_EPSILON
 #define RNMX (1.0-EPS)
 */

double ran2(void)
{
	long j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	
	if (*idum <= 0) 
	{ 
		if (-(*idum) < 1) *idum=1; /*Be sure to prevent idum = 0.*/
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+32;j>=0;j--) 
		{
			/*Load the shuffle table (after 32 warm-ups).*/
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1; /*Start here when not initializing.*/
	/*Compute idum=(IA1*idum) % IM1 without*/
	/*overflows by Schrage's method. */
	*idum=IA1*(*idum-k*IQ1)-k*IR1; 
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	/*Compute idum2=(IA2*idum) % IM2 likewise.*/
	idum2=IA2*(idum2-k*IQ2)-k*IR2; 
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV2; /*Will be in the range 0..NTAB-1.*/
	iy=iv[j]-idum2; 
	/*Here idum is shuffled, idum and idum2 are
	 combined to generate output. */
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	/*the output will be (0, 1]*/
	return  (double)AM2*iy;
	/*	
	 Here the output will be (0, 1)
	 if ((temp=AM2*iy)>RNMX) return RNMX
	 else return temp;
	 */
	
}
#undef IM1
#undef IM2
#undef AM2
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV2
/*--------------------------------------------------------------------------------*/
/*Function returns the poisson deviate given a mean value xm. Original function from 
 the Numerical Recipes in C IIed. modified to use global idum 
 */
double poidev(double xm)
{
	static double sq, alxm, g, oldm=(-1.0);
	double em, t, y;
	if (xm < 12.0) 
	{
		if (xm != oldm) 
		{
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do 
		{
        	++em;
        	t *= ran2();
		} while (t > g);
	} 
	else 
	{
		if (xm != oldm) 
		{
        	oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do 
		{
        	do 
        	{
				y=tan(PI*ran2());
				em=sq*y+xm;
			} while (em < 0.0);
			em=(long)(em);
    		t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran2() > t);
	}
	return em;
}
/*--------------------------------------------------------------------------------*/
/*ln gamma(x) function from the Numerical Recipes in C IIed. */
double gammln(double xx)
{
	double x, tmp, ser;
	static double cof[6]={76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003e-2, -0.536382e-5};
	long j;
	
	x = xx-1.0;
	tmp = x + 5.5;
	tmp -= (x+0.5) * log(tmp);
	ser = 1.0;
	for (j=0; j<=5; j++) 
	{
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}
/*--------------------------------------------------------------------------------*/
/*function to find the time. starts if called with 
 parameter 0. if called with parameter 1 it gets the 
 total time and exits. if called with parameter 2 it will 
 find the time since the last call and keep the counter on*/
long getTime(long m)
{
	static long start=0, end=0;
	static long sttemp=0;
	long diff;
	if(m==0)
	{
		start=time(NULL);
	}
	/*time from start*/
	else if(m==1)
	{
		end=time(NULL);
		diff=difftime(end, start);
		printf("\nTotal Time from start:  %02ld:%02ld:%02ld = %02ld:%02ld\n", diff/3600, (diff%3600)/60, diff%60, diff/60, diff%60);
	}
	/*time since last call*/
	else if(m==2)
	{
		if(end)sttemp=end;
		else sttemp=start;
		end=time(NULL);
		diff=difftime(end, sttemp);
		printf("Last Interval:  %02ld:%02ld:%02ld = %02ld:%02ld\n", diff/3600, (diff%3600)/60, diff%60, diff/60, diff%60);
	}
	/*time from start*/
	else if(m==3)
	{
		end=time(NULL);
		diff=difftime(end, start);
		return diff;
	}
	/*time since last call*/
	else if(m==4)
	{
		if(end)sttemp=end;
		else sttemp=start;
		end=time(NULL);
		diff=difftime(end, sttemp);
		return diff;
	}
	return 0;
}
/*---------------------------------------------------------------------------------------*/



