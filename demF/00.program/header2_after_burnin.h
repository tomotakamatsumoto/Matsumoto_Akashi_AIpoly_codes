/* Period parameters */
#define L 624
#define O 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[L]; /* the array for the state vector  */
static int mti=L+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<L; mti++) {
        mt[mti] =
        (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}
/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (L>key_length ? L : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
        + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=L) { mt[0] = mt[L-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=L-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
        - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=L) { mt[0] = mt[L-1]; i=1; }
    }
    
    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    
    if (mti >= L) { /* generate N words at one time */
        int kk;
        
        if (mti == L+1)   /* if init_genrand() has not been called, */
        init_genrand(5489UL); /* a default initial seed is used */
        
        for (kk=0;kk<L-O;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+O] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<L-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(O-L)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[L-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[L-1] = mt[O-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        
        mti = 0;
    }
    
    y = mt[mti++];
    
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    
    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */



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
#define MAXSEQLEN		300000				/*maximum length of sequence has to be multiple of 3*/
#define MAXSAMPLENUM	10000					/*maximum samplesize*/
#define MAXSEQNUM		10000				/*maximum number of sequences even when fluctuating*/
#define MINSEQNUM		1000					/*minimum number of sequences even under population fluctuation*/
#define MAXGEN			500000				/*maximum of inirungen, prerungen and reprungen*/
#define MAXFIXATIONS	1000000				/*maximum number of fixations permitted in a single initrun or prerun or reprun*/



/*---------------------------------------------------------------------------------------*/
/*structure definitions*/
struct seq_dat_str
{
	char site[MAXSEQLEN + 1];					/*string of bases for each sequence*/
	long baseCount0[3+1];			/*baseCount0[j] is the number of 0s at sitePos where (sitePos+2)%3+1=j*/
	long baseCount1[3+1];           /*baseCount1[j] is the number of 1s at sitePos where (sitePos+2)%3+1=j*/
    long baseCount2[3+1];           /*baseCount2[j] is the number of 2s at sitePos where (sitePos+2)%3+1=j*/
    long baseCount3[3+1];           /*baseCount3[j] is the number of 3s at sitePos where (sitePos+2)%3+1=j*/
    long tot_1_ct;
	long tot_1_ct_silent;
    long tot_1_ct_replacement;
}
*seq_dat, *temp_seq_dat, *prev_seq_dat;

struct siteData_str
{	
	char ancSite;							/*ancestral base at current site*/
	long freqs0;					/*frequency of 0s in the population at current site*/
	long freqs1;					/*frequency of 1s in the population at current site*/
    long freqs2;					/*frequency of 2s in the population at current site*/
    long freqs3;					/*frequency of 3s in the population at current site*/
	long fixNum01;					/*number of fixations from 0 to 1 at current site in the population during current replicate */
    long fixNum02;					/*number of fixations from 0 to 2 at current site in the population during current replicate */
    long fixNum03;					/*number of fixations from 0 to 3 at current site in the population during current replicate */
    long fixNum10;					/*number of fixations from 1 to 0 at current site in the population during current replicate */
    long fixNum12;					/*number of fixations from 1 to 2 at current site in the population during current replicate */
    long fixNum13;					/*number of fixations from 1 to 3 at current site in the population during current replicate */
    long fixNum20;					/*number of fixations from 2 to 0 at current site in the population during current replicate */
    long fixNum21;					/*number of fixations from 2 to 1 at current site in the population during current replicate */
    long fixNum23;					/*number of fixations from 2 to 3 at current site in the population during current replicate */
    long fixNum30;					/*number of fixations from 3 to 0 at current site in the population during current replicate */
    long fixNum31;					/*number of fixations from 3 to 1 at current site in the population during current replicate */
    long fixNum32;					/*number of fixations from 3 to 2 at current site in the population during current replicate */
    
	long mutNum01;					/*number of mutations from 0 to 1 at current site in the population during current replicate */
    long mutNum02;					/*number of mutations from 0 to 2 at current site in the population during current replicate */
    long mutNum03;					/*number of mutations from 0 to 3 at current site in the population during current replicate */
    long mutNum10;					/*number of mutations from 1 to 0 at current site in the population during current replicate */
    long mutNum12;					/*number of mutations from 1 to 2 at current site in the population during current replicate */
    long mutNum13;					/*number of mutations from 1 to 3 at current site in the population during current replicate */
    long mutNum20;					/*number of mutations from 2 to 0 at current site in the population during current replicate */
    long mutNum21;					/*number of mutations from 2 to 1 at current site in the population during current replicate */
    long mutNum23;					/*number of mutations from 2 to 3 at current site in the population during current replicate */
    long mutNum30;					/*number of mutations from 3 to 0 at current site in the population during current replicate */
    long mutNum31;					/*number of mutations from 3 to 1 at current site in the population during current replicate */
    long mutNum32;					/*number of mutations from 3 to 3 at current site in the population during current replicate */

	long recNum;					/*number of crossovers at current site during current replicate */
}
*siteData;

double *w_dat_static_silent;
double *w_dat_static_replacement;
double *w_dat_static_silent_b;
double *w_dat_static_replacement_b;

double init_MCU;
clock_t start_time, end_time;

/*flags read from the control file*/
//long noMultiHit;								/*1=no multiple hits 0=permit multiple hits*/

/*parameters read from the control file*/
//long repNum;						/*Total number of replicates for the simulation*/

//long initRunGen;					/*number of repRunGen for the initial run */
//long preRunGen;					/*number of repRunGen for the prerun before each reprun*/
//long repRunGen;					/*number of generations for the replicate run*/	

//long initSeqNum;					/*Initial number of sequences*/
//long seqLen;						/*Length of sequences as number of bases*/
//long sampleNum;					/*Number of sequences to be sampled to collect data*/

//double c;									/*per site per generation recombination rate*/
//double u10;									/*the 1->0 mutation rate common for all positions*/
//double u01;									/*the 0->1 mutation rate common for all positions*/
//double Nes2;								/* 2Nes for MCU model */

char *outFileName, *outFileName2, *outFile_segsites, *outFile_MCU;			/*name of output file*/
char *outFileExtn;			/*outputfilename extension*/
//char *ctlFileName;			/*name of control file*/
//char version[11];							/*version string will be appended to filename before extension*/
//char fsimVersion[11];						/*version of fsim. the above is the version of the simulation*/
//char *comments;				/*comments - will be output into the outputfile*/
//char *machineName;			/*machine on which the program is running*/
FILE *fpOut, *fpOut2, *fpOut_segsites, *fpOut_MCU;								/*output file pointer*/

/*variables used for the simulation - data being processed*/
double *sArray;					/*array of s values for the sequence*/

long *preRunPopFluct;		/*array to store the population sizes during each generation in the prerun*/
long *popSizeArray;		/*array used to store the population sizes for each gen for a given cycle*/
// long curPopSize;					/*The current population size at any time*/
long curGenNum;					/*the current generation number*/
long curRepNum;					/*the current replicate*/
long *sampleSeqs;	/*indices of sequences that are to be sampled for the current replicate*/		

/* data for doMutations() */
long *hitSites10;	/*indices of 1 sites to be mutated to 0*/
long *hitSites20;	/*indices of 2 sites to be mutated to 0*/
long *hitSites30;	/*indices of 3 sites to be mutated to 0*/
long *hitSites01;	/*indices of 0 sites to be mutated to 1*/
long *hitSites02;	/*indices of 0 sites to be mutated to 2*/
long *hitSites03;	/*indices of 0 sites to be mutated to 3*/
long *hitSites12;	/*indices of 1 sites to be mutated to 2*/
long *hitSites13;	/*indices of 1 sites to be mutated to 3*/
long *hitSites21;	/*indices of 2 sites to be mutated to 1*/
long *hitSites31;	/*indices of 3 sites to be mutated to 1*/
long *hitSites23;	/*indices of 2 sites to be mutated to 3*/
long *hitSites32;	/*indices of 3 sites to be mutated to 2*/


/* data for doRecombination() */
long *seqIndex;			/*Array to store the indices to pair sequences randomly*/
long *crossOverSite;	/*array to store the positions of crossovers*/

/* data for getNextGenSeqs() and getNextGenSeqs_MCU() */
double *wData;				/*array that stores the weight for each sequence*/
double *freqs;				/*expected frequencies of current individuals*/
long *newFreqs;			/*obtained frequencies of offsprings for current generation*/

/*arrays of flags for processing data*/
long *curSeqHitSite; 			/*1=site already hit in current mutagenesis of a gene, 0=not*/
long *curGenHitSite;			/*1=site already hit in current generation (across genes) 0=not*/
long *polySite;					/*1=site is segregating in population, 0=not*/		
long *crossHitSite;				/*1=site already selected as crossover site in a given generation, 0=not*/

/*variables for data collection from each replicate. */
long sampleFix01[3+1];				/*Fixations of 0->1 due to sampling (fixed in sample, but not in pop) */
long sampleFix02[3+1];				/*Fixations of 0->2 due to sampling*/
long sampleFix03[3+1];				/*Fixations of 0->3 due to sampling*/
long sampleFix10[3+1];				/*Fixations of 1->0 due to sampling*/
long sampleFix12[3+1];				/*Fixations of 1->2 due to sampling*/
long sampleFix13[3+1];				/*Fixations of 1->3 due to sampling*/
long sampleFix20[3+1];				/*Fixations of 2->0 due to sampling*/
long sampleFix21[3+1];				/*Fixations of 2->1 due to sampling*/
long sampleFix23[3+1];				/*Fixations of 2->3 due to sampling*/
long sampleFix30[3+1];				/*Fixations of 3->0 due to sampling*/
long sampleFix31[3+1];				/*Fixations of 3->1 due to sampling*/
long sampleFix32[3+1];				/*Fixations of 3->2 due to sampling*/

long sampleFixNum;					/*Total number of fixations due to sampling */
long *sampleFixPos;					/*Position of fixed (due to sampling) sites in the sample*/

long sampleSeg01[3+1];				/*segregating sites in the sample with ancSite 0*/
long sampleSeg02[3+1];
long sampleSeg03[3+1];
long sampleSeg10[3+1];				/*segregating sites in the sample with ancSite 1*/
long sampleSeg12[3+1];
long sampleSeg13[3+1];
long sampleSeg20[3+1];				/*segregating sites in the sample with ancSite 2*/
long sampleSeg21[3+1];
long sampleSeg23[3+1];
long sampleSeg30[3+1];				/*segregating sites in the sample with ancSite 3*/
long sampleSeg31[3+1];
long sampleSeg32[3+1];

long sampleSegNum;					/*Total number of segregating sites in the sample*/
long *sampleSegPos;					/*Position of segregating sites in the sample*/

long popRecNum;					/*Total number of crossovers in the population across all sites */

/*for every fixation in the population, store the derived state, the position in the sequence, and the fitness effect */
long popFixNum;					/*number of mutations fixed in the replicate - used to access the following arrays*/
char *fixSite;				/*the site that got fixed ie. if 0 was fixed then 0 else 1*/
long *fixSitePos;	/*Positions of fixed sites in current replicate*/
double *fixSiteSval;				/*S values of sites fixed in current replicate*/

long repTime;						/*time in seconds for the current replicate*/

long genMultiHitCount;				/*Total of multiple hits in a generation across all generations for current rep*/
long segMultiHitCount;				/*Total of multiple hits at segregating sites across all generations*/

/*function declarations*/
long checkParameters(struct userPref_str *userPrefp);
long writeFsimInfo(struct userPref_str *userPrefp);
long initData(struct userPref_str *userPrefp);
long doFsim(struct userPref_str *userPrefp, long genNum, long flagReset, long repType);
long getPopScenario(struct userPref_str *userPrefp, long repType, long genNum);
long resetData(struct userPref_str *userPrefp);
long doMutations(struct userPref_str *userPrefp, long curGenNum);
long getHitSites(struct userPref_str *userPrefp, char origSite, long baseCount, long curSeq, long codPos, long* hitSites, long numHits);
long getPosNthChar(struct userPref_str *userPrefp, char *site, long codPos, long hitCod, char origSite);
long getPosNthChar_v2(struct userPref_str *userPrefp, char *site, long curSeq, long codPos, long hitCod, char origSite);
long doRecombination(struct userPref_str *userPrefp, long curGenNum);
long getNextGenSeqs(struct userPref_str *userPrefp, long curGenNum);

long updateCounts(struct userPref_str *userPrefp, long curGenNum);
double sampleSval(struct userPref_str *userPrefp, long sitePos, char ancSite);
long writeRepData(struct userPref_str *userPrefp, long curGenNum);
long getRepData(struct userPref_str *userPrefp);
long getSample(struct userPref_str *userPrefp, long curGenNum);
long endSim(struct userPref_str *userPrefp);
long checkBaseCounts(struct userPref_str *userPrefp, long pos, long curGenNum);

void CheckAncSite(struct userPref_str *userPrefp);
void get_static_w_dat(struct userPref_str *userPrefp);
long getNextGenSeqs_MCU(struct userPref_str *userPrefp, long curGenNum);

long getTime_print(struct userPref_str *userPrefp, long m);

/*--------------------------------------------------------------------------------*/
/*function that initializes all the variables. whenever this function is called
the whole data will be reset. used only at the begining of the whole simulation (not for each replicate) 
some of the variables are reset between prerun and replicate */
long initData(struct userPref_str *userPrefp)
{
    init_genrand(time(NULL));
	long debug = 0, ct;
	long numOnes_s, numOnes_r;					/*the number of ones in a sequence*/
	long numZeros_s, numZeros_r;					/*the number of zeros in a sequence*/
	long curSite;					/*refers to the current site*/
	long curSeq;					/*refers to the current sequence*/
	long randSite;					/*site that is randomly picked to be filled with a 1*/
    double randvar;					/*random variable used for binomial sampling*/
	double mcu_s, mcu_r, u_v, initNes;						/*mcu for the sequence*/
	long curCodPos=0;							/*current codon position, 1, 2=rep, 3=sil*/	
//	long *site;						/*temporary array to initialize the sequence*/
	
	printf("Initializing Data\n");
//	site = (long *) memAlloc(MAXSEQLEN + 1, sizeof(long), "");
	
	/* memory allocations for global variables */
	outFileName			= (char *) memAlloc (MAXLINECHAR + 1, sizeof(char), "outFileName");
	outFileName2		= (char *) memAlloc (MAXLINECHAR + 1, sizeof(char), "outFileName2");
	outFile_segsites	= (char *) memAlloc (MAXLINECHAR + 1, sizeof(char), "outFile_segsites");
	outFile_MCU			= (char *) memAlloc (MAXLINECHAR + 1, sizeof(char), "outFile_MCU");
	
	seq_dat			= (struct seq_dat_str *) memAlloc (MAXSEQNUM + 1, sizeof(struct seq_dat_str), "seq array");
    prev_seq_dat	= (struct seq_dat_str *) memAlloc (MAXSEQNUM + 1, sizeof(struct seq_dat_str), "seq array");
	temp_seq_dat	= (struct seq_dat_str *) memAlloc (MAXSEQNUM + 1, sizeof(struct seq_dat_str), "new seq array");
// 	for (curSeq = 0; curSeq <= userPrefp->initSeqNum; curSeq++)
// 	{
// 		seq_dat[curSeq].site	= (char *) memAlloc (userPrefp->seqLen + 1, sizeof(char), "");
// 		temp_seq_dat[curSeq].site	= (char *) memAlloc (userPrefp->seqLen + 1, sizeof(char), "");
// 	}
	siteData = (struct siteData_str *) memAlloc (userPrefp->seqLen + 1, sizeof(struct siteData_str), "siteData"); 

	sArray = (double *) memAlloc (userPrefp->seqLen + 1, sizeof(double), "sArray");

	preRunPopFluct 		= (long *) memAlloc (MAXGEN + 1, sizeof(long), "preRunPopFluct");
	popSizeArray 		= (long *) memAlloc (MAXGEN + 1, sizeof(long), "popSizeArray");
	sampleSeqs 			= (long *) memAlloc (MAXSAMPLENUM + 1, sizeof(long), "popSizeArray");
	
	curSeqHitSite 		= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), ""); 
	curGenHitSite 		= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "");
	polySite 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "");
	crossHitSite 		= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "");
	
	sampleFixPos		= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), ""); 
	sampleSegPos		= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), ""); 

	fixSite				= (char *) memAlloc (MAXFIXATIONS + 1, sizeof(char), "fixSite");
	fixSitePos			= (long *) memAlloc (MAXFIXATIONS + 1, sizeof(long), "fixSitePos");
	fixSiteSval  		= (double *) memAlloc (MAXFIXATIONS + 1, sizeof(double), "fixSiteSval");

	/* memory allocations for doMutations() */
	hitSites01 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "hitSites01");
    hitSites02 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "hitSites02");
    hitSites03 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "hitSites03");
    hitSites10 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "hitSites10");
    hitSites12 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "hitSites12");
    hitSites13 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "hitSites13");
    hitSites20 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "hitSites20");
    hitSites21 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "hitSites21");
    hitSites23 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "hitSites23");
    hitSites30 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "hitSites30");
    hitSites31 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "hitSites31");
    hitSites32 			= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "hitSites32");


	/* memory allocations for doRecombination() */
	seqIndex 			= (long *) memAlloc (userPrefp->prerunSeqNum + 1, sizeof(long), "seqIndex");
	crossOverSite 		= (long *) memAlloc (userPrefp->seqLen + 1, sizeof(long), "crossOverSite");
	
	/* memory allocations for getNextGenSeqs() */
    /* 141205 matsumoto changed the size initSeqNum -> MAXSEQNUM */
	wData 				= (double *) memAlloc (MAXSEQNUM + 1, sizeof(double), "wData");
	freqs 				= (double *) memAlloc (MAXSEQNUM + 1, sizeof(double), "freqs");
	newFreqs 			= (long *) memAlloc (MAXSEQNUM + 1, sizeof(long), "newFreqs");


	

	/* initialize output files */
	sprintf(outFileName,  "%s%s%s_1_%s.%s", userPrefp->rootOutputFolder, userPrefp->dirSep, userPrefp->outFileName, userPrefp->version, userPrefp->outFileExtn);
	sprintf(outFileName2, "%s%s%s_2_%s.%s", userPrefp->rootOutputFolder, userPrefp->dirSep, userPrefp->outFileName, userPrefp->version, userPrefp->outFileExtn);
	sprintf(outFile_segsites, "%s%s%s_segsites_%s.%s", userPrefp->rootOutputFolder, userPrefp->dirSep, userPrefp->outFileName, userPrefp->version, userPrefp->outFileExtn);
	sprintf(outFile_MCU, "%s%s%s_MCU_%s.%s", userPrefp->rootOutputFolder, userPrefp->dirSep, userPrefp->outFileName, userPrefp->version, userPrefp->outFileExtn);
	
	fpOut			= fileOpen(outFileName, "wb");
	fpOut2			= fileOpen(outFileName2, "wb");
	fpOut_segsites	= fileOpen(outFile_segsites, "wb");
	fpOut_MCU		= fileOpen(outFile_MCU, "wb");
	
	idumVal = userPrefp->idum_init;
	checkParameters(userPrefp);
	
	if (debug == 1)	printf("\tpos 1\n");
	if (debug == 1)	printf("\tpos 2\n");
	
	/* calculate the initial MCU */
    if (userPrefp->fitness_u == 1.000) {
        init_MCU = 0.50;
    }
    else if (userPrefp->fitness_p < 1.000) {
        initNes = userPrefp->prerunSeqNum * (1.000 - (userPrefp->fitness_u));
        u_v = 1.000;
        init_MCU = exp(2.0 * initNes);	/* expected MCU for 2f sites under free recombination */
        init_MCU /= (exp(2.0 * initNes) + u_v);
    }
//	mcu = INITMCU;								/* this should not be initialized to 0.50 for MCU model !!!!!!!!!!!!!!!!!!!!!! */
	//mcu_s = init_MCU;
	//numOnes_s = mcu_s * (userPrefp->seqLen/2);
	//numZeros_s = (userPrefp->seqLen/2) - numOnes_s;
    //mcu_r = userPrefp->u01_a / (userPrefp->u10_a + userPrefp->u01_a);
    //numOnes_r = mcu_r * (userPrefp->seqLen/2);
    //numZeros_r = (userPrefp->seqLen/2) - numOnes_r;
	
	/*iniitializes the temporary array and fills the first sequence with all zeros also initializes the siteData structure*/
    
    /*150624 matsumoto: set the initial population (all individuals has a sequence of codon TTT)*/
    
    for (curSeq=1; curSeq<= userPrefp->initSeqNum; curSeq++) {
        for(curSite = 1; curSite <= userPrefp->seqLen; curSite++) {
            seq_dat[curSeq].site[curSite] = userPrefp->init_seq_dat[curSeq][curSite-1];
        }
    }
    for (curSeq=1; curSeq<= userPrefp->initSeqNum; curSeq++)
    {
        seq_dat[curSeq].baseCount0[1] = seq_dat[curSeq].baseCount0[2] = seq_dat[curSeq].baseCount0[3] = 0;
        seq_dat[curSeq].baseCount1[1] = seq_dat[curSeq].baseCount1[2] = seq_dat[curSeq].baseCount1[3] = 0;
        seq_dat[curSeq].baseCount2[1] = seq_dat[curSeq].baseCount2[2] = seq_dat[curSeq].baseCount2[3] = 0;
        seq_dat[curSeq].baseCount3[1] = seq_dat[curSeq].baseCount3[2] = seq_dat[curSeq].baseCount3[3] = 0;
        
        for(curSite = 1; curSite <= userPrefp->seqLen; curSite++) {
            if (seq_dat[curSeq].site[curSite] == '1') {
                seq_dat[curSeq].baseCount1[sitePosToCodPos(curSite)]++;
                siteData[curSite].freqs1 = siteData[curSite].freqs1 + 1;
                siteData[curSite].freqs0 = siteData[curSite].freqs0;
                siteData[curSite].freqs2 = siteData[curSite].freqs2;
                siteData[curSite].freqs3 = siteData[curSite].freqs3;
            }
            else if (seq_dat[curSeq].site[curSite] == '2') {
                seq_dat[curSeq].baseCount2[sitePosToCodPos(curSite)]++;
                siteData[curSite].freqs2 = siteData[curSite].freqs2 + 1;
                siteData[curSite].freqs0 = siteData[curSite].freqs0;
                siteData[curSite].freqs1 = siteData[curSite].freqs1;
                siteData[curSite].freqs3 = siteData[curSite].freqs3;
            }
            else if (seq_dat[curSeq].site[curSite] == '3') {
                seq_dat[curSeq].baseCount3[sitePosToCodPos(curSite)]++;
                siteData[curSite].freqs3 = siteData[curSite].freqs3 + 1;
                siteData[curSite].freqs0 = siteData[curSite].freqs0;
                siteData[curSite].freqs1 = siteData[curSite].freqs1;
                siteData[curSite].freqs2 = siteData[curSite].freqs2;
            }
            else if (seq_dat[curSeq].site[curSite] == '0') {
                seq_dat[curSeq].baseCount0[sitePosToCodPos(curSite)]++;
                siteData[curSite].freqs0 = siteData[curSite].freqs0 + 1;
                siteData[curSite].freqs2 = siteData[curSite].freqs2;
                siteData[curSite].freqs1 = siteData[curSite].freqs1;
                siteData[curSite].freqs3 = siteData[curSite].freqs3;
            }
        }
    }

    
    for(curSite = 1; curSite <= userPrefp->seqLen; curSite++)
    {
        if(userPrefp->anc_init_seq_dat[curSite-1] == '1') {
           siteData[curSite].ancSite = '1';
        }
        else if(userPrefp->anc_init_seq_dat[curSite-1] == '2') {
            siteData[curSite].ancSite = '2';
        }
        else if(userPrefp->anc_init_seq_dat[curSite-1] == '3') {
            siteData[curSite].ancSite = '3';
        }
        else if(userPrefp->anc_init_seq_dat[curSite-1] == '0') {
            siteData[curSite].ancSite = '0';
        }
        
        /*Reset Counts in siteData structure*/
        siteData[curSite].fixNum01 = 0;
        siteData[curSite].fixNum10 = 0;
        siteData[curSite].fixNum02 = 0;
        siteData[curSite].fixNum20 = 0;
        siteData[curSite].fixNum03 = 0;
        siteData[curSite].fixNum30 = 0;
        siteData[curSite].fixNum12 = 0;
        siteData[curSite].fixNum21 = 0;
        siteData[curSite].fixNum13 = 0;
        siteData[curSite].fixNum31 = 0;
        siteData[curSite].fixNum23 = 0;
        siteData[curSite].fixNum32 = 0;
        
        siteData[curSite].mutNum01 = 0;
        siteData[curSite].mutNum10 = 0;
        siteData[curSite].mutNum02 = 0;
        siteData[curSite].mutNum20 = 0;
        siteData[curSite].mutNum03 = 0;
        siteData[curSite].mutNum30 = 0;
        siteData[curSite].mutNum12 = 0;
        siteData[curSite].mutNum21 = 0;
        siteData[curSite].mutNum13 = 0;
        siteData[curSite].mutNum31 = 0;
        siteData[curSite].mutNum23 = 0;
        siteData[curSite].mutNum32 = 0;
        
        siteData[curSite].recNum = 0;
    }

    if (debug == 1)	printf("\tpos 4\n");
    
    //	printf("preRunGen 12: %0ld\n", preRunGen);
    /*check if everything was ok*/
    //if (seq_dat[curSeq].baseCount0[1] + seq_dat[curSeq].baseCount0[2] + seq_dat[curSeq].baseCount0[3] != numZeros)
    //errorOut(("Invalid 0 counts(!=%lu) while initializing sequence", numZeros));
    //if (seq_dat[curSeq].baseCount1[1] + seq_dat[curSeq].baseCount1[2] + seq_dat[curSeq].baseCount1[3] != numOnes)
    //errorOut(("Invalid 1 counts(!=%lu) while initializing sequence", numOnes));
    
	
//	printf("preRunGen 12a: %0ld\n", preRunGen);
	/*copy seq_dat[1] to all sequences*/
    /*for (curSeq=2; curSeq<= userPrefp->initSeqNum; curSeq++)
    {

        for (curSite=1; curSite <= userPrefp->seqLen; curSite++) {
            seq_dat[curSeq].site[curSite] = seq_dat[1].site[curSite];
        }

        //seq_dat[curSeq].baseCount0[curSite] = seq_dat[1].baseCount0[curSite];
        //seq_dat[curSeq].baseCount1[curSite] = seq_dat[1].baseCount1[curSite];
        //seq_dat[curSeq].baseCount2[curSite] = seq_dat[1].baseCount2[curSite];
        //seq_dat[curSeq].baseCount3[curSite] = seq_dat[1].baseCount3[curSite];
        
        for (curCodPos = 1; curCodPos <= 3; curCodPos++)
        {
            seq_dat[curSeq].baseCount0[curCodPos] = seq_dat[1].baseCount0[curCodPos];
            seq_dat[curSeq].baseCount1[curCodPos] = seq_dat[1].baseCount1[curCodPos];
            seq_dat[curSeq].baseCount2[curCodPos] = seq_dat[1].baseCount2[curCodPos];
            seq_dat[curSeq].baseCount3[curCodPos] = seq_dat[1].baseCount3[curCodPos];
        }

        
        seq_dat[curSeq].tot_1_ct = seq_dat[1].tot_1_ct;
    }*/
    if (debug == 1)
    {
        for (curSeq=1; curSeq<=userPrefp->initSeqNum; curSeq++)
        {
            printf("seq: %5ld\n", curSeq);
            for (curSite = 1; curSite <= userPrefp->seqLen; curSite++)
            printf("%c", seq_dat[curSeq].site[curSite]);
            printf("\n");
        }
        printf("\n");
    }
    //for (curSite=1; curSite <= userPrefp->seqLen; curSite++) {
    //    printf ("%c", seq_dat[10].site[curSite]);
    //}
    //printf("\n");
	
//	printf("preRunGen 12b: %0ld\n", preRunGen);
	/* Initializes the sArray */
	for (curSite=1; curSite<=userPrefp->seqLen; curSite++)
	{
		/* sArray[curSite]=sampleSval(curSite, seq_dat[1].site[curSite]); */
		sArray[curSite]=0.0;
	}
//	printf("preRunGen 12c: %0ld\n", preRunGen);
	for (curSite=1; curSite<=MAXGEN; curSite++)
	{
		popSizeArray[curSite] = userPrefp->prerunSeqNum;
	}
//	curPopSize = userPrefp->initSeqNum;

//	printf("preRunGen 12d: %0ld\n", preRunGen);
	for(curSeq=1; curSeq <= MAXSAMPLENUM; curSeq++)
		sampleSeqs[curSeq]=0;
	
	if (debug == 1)	printf("\tpos 5\n");
//	printf("preRunGen 13: %0ld\n", preRunGen);

	/*Initializes the flags*/
	for (curSite=1; curSite<=userPrefp->seqLen; curSite++)
	{
		curSeqHitSite[curSite] = 0;
		curGenHitSite[curSite] = 0;
        
        /* 151221 matsumoto: check polymorphic site in the initial populaitn */
        polySite[curSite] = 0;
        for (curSeq=2; curSeq<= userPrefp->initSeqNum; curSeq++) {
            if (seq_dat[curSeq].site[curSite] != seq_dat[1].site[curSite]) {
                polySite[curSite] = 1;
            }
        }
		crossHitSite[curSite] = 0;
	}
	
	/*Initializes the counts*/
    sampleFix01[1]=sampleFix01[2]=sampleFix01[3]=0;
    sampleFix02[1]=sampleFix02[2]=sampleFix02[3]=0;
    sampleFix03[1]=sampleFix03[2]=sampleFix03[3]=0;
    sampleFix10[1]=sampleFix10[2]=sampleFix10[3]=0;
    sampleFix12[1]=sampleFix12[2]=sampleFix12[3]=0;
    sampleFix13[1]=sampleFix13[2]=sampleFix13[3]=0;
    sampleFix20[1]=sampleFix20[2]=sampleFix20[3]=0;
    sampleFix21[1]=sampleFix21[2]=sampleFix21[3]=0;
    sampleFix23[1]=sampleFix23[2]=sampleFix23[3]=0;
    sampleFix30[1]=sampleFix30[2]=sampleFix30[3]=0;
    sampleFix31[1]=sampleFix31[2]=sampleFix31[3]=0;
    sampleFix32[1]=sampleFix32[2]=sampleFix32[3]=0;
    sampleFixNum=0;
    
    sampleSeg01[1]=sampleSeg01[2]=sampleSeg01[3]=0;
    sampleSeg02[1]=sampleSeg02[2]=sampleSeg02[3]=0;
    sampleSeg03[1]=sampleSeg03[2]=sampleSeg03[3]=0;
    sampleSeg10[1]=sampleSeg10[2]=sampleSeg10[3]=0;
    sampleSeg12[1]=sampleSeg12[2]=sampleSeg12[3]=0;
    sampleSeg13[1]=sampleSeg13[2]=sampleSeg13[3]=0;
    sampleSeg20[1]=sampleSeg20[2]=sampleSeg20[3]=0;
    sampleSeg21[1]=sampleSeg21[2]=sampleSeg21[3]=0;
    sampleSeg23[1]=sampleSeg23[2]=sampleSeg23[3]=0;
    sampleSeg30[1]=sampleSeg30[2]=sampleSeg30[3]=0;
    sampleSeg31[1]=sampleSeg31[2]=sampleSeg31[3]=0;
    sampleSeg32[1]=sampleSeg32[2]=sampleSeg32[3]=0;
    sampleSegNum=0;
    
	for (curSite=1;curSite<=userPrefp->seqLen;curSite++)
	{
		sampleFixPos[curSite]=0;
		sampleSegPos[curSite]=0;
	}
	popRecNum=0;
	popFixNum=0;
	for(curSite=1; curSite<=MAXFIXATIONS; curSite++)
	{
		fixSite[curSite]     = '-';
		fixSitePos[curSite]  = 0;
		fixSiteSval[curSite] = 0;
	}
	repTime=0;
	genMultiHitCount=0;
	segMultiHitCount=0;
	if (debug == 1)	printf("\tpos 9\n");
	
//	printf("preRunGen 7: %0ld\n", preRunGen);
//	free(site);
	return 0;
}
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/*Checks the input parameters*/
long checkParameters(struct userPref_str *userPrefp)
{	
	long maxGen;					/*maximum length of any interval*/
	long expFixations;				/*expected fixations during the max interval*/
	long fourZigmaRange;			/*four zigma range for the number of fixations*/
	double maxU;							/*maximum of the two mutation rates*/
//	long i, j, imax;							/*temporary ints for looping*/
	
	/*Program mode*/
	#ifdef DEBUGMODE
		printf("Program Running in DEBUGMODE mode.\n"
				"Could be slower and there could be more outputs than needed\n"
				"Undefine DEBUGMODE to change");
	#endif
	#ifdef VERBOSE
		printf("Program Running in VERBOSE mode.\n"
				"Could be slower and there could be more outputs than needed\n"
				"Undefine VERBOSE to change");
	#endif
	
	/*checking filename version number and output files*/
	/*check if version inside ctl file matches version in ctl file name*/
//	j=strlen(userPrefp->version);
//	imax=strlen(ctlFileName)-j+1;
//	i=0;
//	while(i<imax)
//	{
//		if(!strncmp(version, &ctlFileName[i], j))
//			break;
//		i++;
//	}
//	if(i==imax)
//		errorOut(("version number in control file name doesnt match version number inside file(%s)\n", version));
	
	/*check if output file already exists*/
// 	fpOut=fopen(outFileName, "rb");
// 	if (fpOut)
// 		errorOut(("output file '%s' exists. delete it or rename it", outFileName));
// 	else
// 		fclose(fpOut);
// 	
// 	fpOut2=fopen(outFileName2, "rb");
// 	if(fpOut2)
// 		errorOut(("output file '%s' exists. delete it or rename it", outFileName2));
// 	else
// 		fclose(fpOut2);
// 	
// 	fpOut_segsites = fopen(outFile_segsites, "rb");
// 	if(fpOut_segsites)
// 		errorOut(("output file '%s' exists. delete it or rename it", outFile_segsites));
// 	else
// 		fclose(fpOut_segsites);
// 	
// 	fpOut_MCU= fopen(outFile_MCU, "rb");
// 	if(fpOut_MCU)
// 		errorOut(("output file '%s' exists. delete it or rename it", outFile_MCU));
// 	else
// 		fclose(fpOut_MCU);

	/*population size*/
	if(userPrefp->initSeqNum > MAXSEQNUM || userPrefp->initSeqNum < MINSEQNUM)
		errorOut(("initSeqNum(%lu) should be between MAXSEQNUM(%ld) and MINSEQNUM(%ld)", userPrefp->initSeqNum, (long) MAXSEQNUM, (long) MINSEQNUM));
	if(userPrefp->initSeqNum % 2)
		errorOut(("initSeqNum(%lu) has to be a multiple of 2", userPrefp->initSeqNum));
	if(MAXSEQNUM % 2)
		errorOut(("MAXSEQNUM(%lu) has to be a multiple of 2", (long) MAXSEQNUM));
	if(MINSEQNUM % 2)
		errorOut(("MAXSEQNUM(%lu) has to be a multiple of 2", (long) MINSEQNUM));
	
	/*sample size*/
	
	if(userPrefp->sampleNum>MAXSAMPLENUM)
		errorOut(("sampleSize(%lu) cannot be greater than MAXSAMPLENUM(%lu)", userPrefp->sampleNum, (long) MAXSAMPLENUM));
	
	/*generation sizes*/
	if (userPrefp->initRunGen > MAXGEN || userPrefp->preRunGen > MAXGEN || userPrefp->repRunGen > MAXGEN)
		errorOut(("initRunGen, preRunGen and repRunGen should be smaller than MAXGEN(%lu)", (long) MAXGEN));
	
	/*sequence lengths*/
	if (userPrefp->seqLen > MAXSEQLEN)
		errorOut(("seqLen(%lu) cannot be greater than MAXSEQLEN(%lu)", userPrefp->seqLen, (long) MAXSEQLEN));
	//if (userPrefp->seqLen % 3)
	//	errorOut(("seqLen(%lu) has to be a multiple of 3", userPrefp->seqLen));
	if (MAXSEQLEN % 3)
		errorOut(("MAXSEQLEN(%lu) has to be a multiple of 3", (long) MAXSEQLEN));
	
	/*fixations per interval*/
	maxU = (userPrefp->u_a > userPrefp->u_b) ? userPrefp->u_a: userPrefp->u_b;
	maxGen = (userPrefp->initRunGen > userPrefp->preRunGen) ? userPrefp->initRunGen : userPrefp->preRunGen;
	maxGen = (userPrefp->repRunGen > maxGen) ? userPrefp->repRunGen : maxGen;
	expFixations = maxGen * maxU * userPrefp->seqLen;
	fourZigmaRange = expFixations + 4 * sqrt(expFixations);
	if(MAXFIXATIONS < fourZigmaRange)
		errorOut(("MAXFIXATIONS(%lu) less than four Zigma Range(%lu)", (long) MAXFIXATIONS, fourZigmaRange));
	
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*writes the simulation information as the main header in the output file*/
long writeFsimInfo(struct userPref_str *userPrefp)
{
//	printf("preRunGen 6: %0ld\n", preRunGen);

	fprintf(fpOut, "fileName           = %s\n", 		outFileName		);
	fprintf(fpOut, "version            = %s\n", 		userPrefp->version			);
	fprintf(fpOut, "fsimVersion        = %s\n", 		userPrefp->fsimVersion		);
	fprintf(fpOut, "machineName        = %s\n", 		userPrefp->machineName		);
	fprintf(fpOut, "MAXSEQNUM          = %d\n", 		MAXSEQNUM		);
	fprintf(fpOut, "MINSEQNUM          = %d\n", 		MINSEQNUM		);
	fprintf(fpOut, "POP_MODEL          = %s\n", 		POP_MODEL		);
    fprintf(fpOut, "prerunSeqNum         = %ld\n", 		userPrefp->prerunSeqNum		);
	fprintf(fpOut, "initSeqNum         = %ld\n", 		userPrefp->initSeqNum		);
    fprintf(fpOut, "SeqNum_b         = %ld\n", 		userPrefp->SeqNum_b		);
	fprintf(fpOut, "sampleNum          = %ld\n", 		userPrefp->sampleNum		);
	fprintf(fpOut, "seqLen             = %ld\n", 		userPrefp->seqLen			);
	fprintf(fpOut, "initRunGen         = %ld\n", 		userPrefp->initRunGen		);
	fprintf(fpOut, "preRunGen          = %ld\n", 		userPrefp->preRunGen		);
	fprintf(fpOut, "repRunGen          = %ld\n", 		userPrefp->repRunGen		);
	fprintf(fpOut, "repNum             = %ld\n", 		userPrefp->repNum			);
	fprintf(fpOut, "u_a                = %.10f\n", 		userPrefp->u_a				);
	fprintf(fpOut, "u_b                = %.10f\n", 		userPrefp->u_b				);
	fprintf(fpOut, "c                  = %.10f\n", 		userPrefp->c				);
    fprintf(fpOut, "c_b                  = %.10f\n", 		userPrefp->c_b				);
	fprintf(fpOut, "noMultiHit         = %ld\n", 		userPrefp->noMultiHit		);
	fprintf(fpOut, "S_MODEL            = %s\n", 		S_MODEL			);
	fprintf(fpOut, "idum               = %ld\n", 		userPrefp->idum_init		);
	fprintf(fpOut, "fsimComments       = %s\n", 		userPrefp->comments		);
	fprintf(fpOut, "//\n");
	fflush(fpOut);
	
	printf("Writing simulation information to output file\n");
	printf("fileName           = %s\n", 			outFileName		);
	printf("version            = %s\n", 			userPrefp->version			);
	printf("fsimVersion        = %s\n", 			userPrefp->fsimVersion		);
	printf("machineName        = %s\n", 			userPrefp->machineName		);
	printf("MAXSEQNUM          = %d\n", 			MAXSEQNUM		);
	printf("MINSEQNUM          = %d\n", 			MINSEQNUM		);
	printf("POP_MODEL          = %s\n", 			POP_MODEL		);
    printf("prerunSeqNum       = %ld\n", 		userPrefp->prerunSeqNum		);
	printf("initSeqNum         = %ld\n", 			userPrefp->initSeqNum		);
	printf("sampleNum          = %ld\n", 			userPrefp->sampleNum		);
	printf("seqLen             = %ld\n", 			userPrefp->seqLen			);
	printf("initRunGen         = %ld\n", 			userPrefp->initRunGen		);
	printf("preRunGen          = %ld\n", 			userPrefp->preRunGen		);
	printf("repRunGen          = %ld\n", 			userPrefp->repRunGen		);
	printf("repNum             = %ld\n", 			userPrefp->repNum			);
	printf("u_a                = %.10f\n", 			userPrefp->u_a				);
	printf("u_b                = %.10f\n", 			userPrefp->u_b				);
	printf("c                  = %.10f\n", 			userPrefp->c				);
	printf("noMultiHit         = %ld\n", 			userPrefp->noMultiHit		);
	printf("S_MODEL            = %s\n", 			S_MODEL			);
	printf("idum               = %ld\n", 			userPrefp->idum_init		);
	printf("fsimComments       = %s\n", 			userPrefp->comments		);
	halt;
	
	return 0;
}
/*--------------------------------------------------------------------------------*/
/* 
	output program run time to file in 20.output folder
	original code from anoop gettime() but this has been altered
	the function is currently only defined for m=1 (the other options are left for future coding)

 declare variables: 
 	clock_t start_time, end_time;
 
 call clock functions:
 	start_time = clock();
	end_time = clock();
	
 then call this function

*/
long getTime_print(struct userPref_str *userPrefp, long m)
{
	FILE *fpoutputfile;
	char outputfile[MAXFILENAME + 1];
	static long start=0, end=0;
	static long sttemp=0;
	long diff, secs_int;
	double secs = 0;
	
	sprintf(outputfile, "%s%s%s", userPrefp->rootOutputFolder, userPrefp->dirSep, "0.time_out");
	fpoutputfile = fileOpen(outputfile, "w");

	if(m==0)
	{
		start=time(NULL);
	}
	/*time from start*/
	else if(m==1)
	{
// 		end=time(NULL);
// 		diff=difftime(end, start);
// 		fprintf(fpoutputfile, "\nTotal Time from start:  %02ld:%02ld:%02ld = %02ld:%02ld\n", diff/3600, (diff%3600)/60, diff%60, diff/60, diff%60);
// 		fprintf(fpoutputfile, "\nTotal sec from start:  %0ld\n", diff);
		secs = (double) (end_time - start_time) / CLOCKS_PER_SEC; 
		secs_int = (long) secs;
		fprintf(fpoutputfile, "\nElapsed seconds = %0.0f\n", secs);
		fprintf(fpoutputfile, "\nTotal Time from start:  %02ld:%02ld:%02ld = %02ld:%02ld\n", secs_int/3600, (secs_int%3600)/60, secs_int%60, secs_int/60, secs_int%60);
		fprintf(fpoutputfile, "\nTotal sec from start:  %0ld\n", secs_int);
	}
	/*time since last call*/
	else if(m==2)
	{
		if(end)sttemp=end;
		else sttemp=start;
		end=time(NULL);
		diff=difftime(end, sttemp);
		fprintf(fpoutputfile, "Last Interval:  %02ld:%02ld:%02ld = %02ld:%02ld\n", diff/3600, (diff%3600)/60, diff%60, diff/60, diff%60);
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
	
	fclose(fpoutputfile);

	return 0;
}
/*--------------------------------------------------------------------------------*/



