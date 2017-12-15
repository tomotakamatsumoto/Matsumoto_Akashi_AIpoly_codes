/*
FsimII
-------
Program created :	05/15/03	Anoop John
Last modified	:	07/21/03	Anoop John

The  second version of fsim created from  scratch to  rectify
problems  and  limitations  with  old fsim. The program  will
have to  be simple,  easy  to  modify  and  easy  to upgrade.
Simplicity  will  be  a  higher priority   than  performance.
The  program  will  read  in  parameters from  a control file
and create an output file at the given  location. Output file
format will be maintained in a similar manner as the previous
fsim  to be able to use the old read program with very little
modifications. The last version of old fsim was fsim_AJ_080.

version.00
Started	:	05/15/03
version.01
Started :	05/27/03
version.02
Started :	06/04/03
version.03
Started :	06/06/03
version.04
Started :	06/16/03

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
			For 1,2,3 codon positions
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
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
#define INITMCU			0.5					/*Initial MCU value used to initialize the sequences                                  CHANGE THIS TO EQUILIBRIUM VALUE FOR MCU MODEL */
#define S_MODEL			"sampleAtEveryMutation"	/*sampling model for s values*/

/*constants*/
#define MAXSEQLEN		1200				/*maximum length of sequence has to be multiple of 3*/
#define MAXSAMPLENUM	25					/*maximum samplesize*/	
#define	MAXSEQNUM		5000ul				/*maximum number of sequences even when fluctuating*/
#define MINSEQNUM		500ul				/*minimum number of sequences even under population fluctuation*/
#define MAXGEN			100000				/*maximum of inirungen,prerungen and reprungen*/
#define MAXFIXATIONS	200000				/*maximum number of fixations permitted in a single initrun or prerun or reprun*/
#define FLOATPRECISION	0.000001			/*precision to be used for calculations*/
#define PI				3.1415926535897932384626433832795
#define MAXLINECHAR		1000				/*maximum characters in a line*/
#define	MAXFILENAME		31					/*maximum length of filename*/
#define IDUMINIT		-123341l			/*value used to initialize the random number generator*/
/*---------------------------------------------------------------------------------------*/
/*structure definitions*/
struct seq_str
{
	char site[MAXSEQLEN+1];					/*string of bases for each sequence*/
	unsigned long baseCount0[3+1];			/*baseCount0[j] is the number of 0s at sitePos where (sitePos+2)%3+1=j*/
	unsigned long baseCount1[3+1];			/*baseCount1[j] is the number of 1s at sitePos where (sitePos+2)%3+1=j*/
};
struct siteData_str
{	
	char ancSite;							/*ancestral base at current site*/
	unsigned long freqs0;					/*frequency of 0s in the population at current site*/
	unsigned long freqs1;					/*frequency of 1s in the population at current site*/
	unsigned long fixNum01;					/*number of fixations from 0 to 1 at current site in the population during current replicate */
	unsigned long fixNum10;					/*number of fixations from 1 to 0 at current site in the population during current replicate */
	unsigned long mutNum01;					/*number of mutations from 0 to 1 at current site in the population during current replicate */
	unsigned long mutNum10;					/*number of mutations from 1 to 0 at current site in the population during current replicate */
	unsigned long recNum;					/*number of crossovers at current site during current replicate */
};
/*---------------------------------------------------------------------------------------*/
/*debug tools*/
#define pf(a) 		printf("%.20f\t",a)
#define pl(a) 		printf("%ld\t",a)
#define ps(a) 		printf("%s\t",a)
#define pa(a) 		printf("%c\t",a)
#define tab			printf("\t")
#define ent			printf("\n")
#define help		printf(" O'er h're")
#define halt		do{printf("\nPress enter to continue!\n");getchar();}while(0)
//#undef halt 
//#define halt		do{printf("\nPress any key to continue!\n");}while(0)
/*---------------------------------------------------------------------------------------*/
/*Error function implemented as macro*/
#define errorOut(a)	do{printf("\nError! ");printf a;printf("\nexiting program.....");getchar();exit(1);}while(0)
/*---------------------------------------------------------------------------------------*/
/*macro function to find codon position given site position*/
#define sitePosToCodPos(a)	((((a)+2)%3)+1)
/*---------------------------------------------------------------------------------------*/
/*global variables*/

/*flags read from the control file*/
int noMultiHit;								/*1=no multiple hits 0=permit multiple hits*/

/*parameters read from the control file*/
unsigned long repNum;						/*Total number of replicates for the simulation*/

unsigned long initRunGen;					/*number of repRunGen for the initial run */
unsigned long preRunGen;					/*number of repRunGen for the prerun before each reprun*/
unsigned long repRunGen;					/*number of generations for the replicate run*/	

unsigned long initSeqNum;					/*Initial number of sequences*/
unsigned long seqLen;						/*Length of sequences as number of bases*/
unsigned long sampleNum;					/*Number of sequences to be sampled to collect data*/

double c;									/*per site per generation recombination rate*/
double u10;									/*the 1->0 mutation rate common for all positions*/
double u01;									/*the 0->1 mutation rate common for all positions*/

long idumVal=IDUMINIT;						/*default ran2 seed. if seed given!=0 then it will be used*/	
long *idum=&idumVal;						
char outFileName[MAXLINECHAR+1];			/*name of output file*/
char outFileExtn[MAXLINECHAR+1];			/*outputfilename extension*/
char ctlFileName[MAXLINECHAR+1];			/*name of control file*/
char version[11];							/*version string will be appended to filename before extension*/
char fsimVersion[11];						/*version of fsim. the above is the version of the simulation*/
char comments[MAXLINECHAR+1];				/*comments - will be output into the outputfile*/
char machineName[MAXLINECHAR+1];			/*machine on which the program is running*/
FILE *fpOut;								/*output file pointer*/

/*variables used for the simulation - data being processed*/
struct seq_str *seq,*newSeq;
struct siteData_str	siteData[MAXSEQLEN+1];	/*array of sitedata structures */
double sArray[MAXSEQLEN+1];					/*array of s values for the sequence*/

unsigned long preRunPopFluct[MAXGEN+1];		/*array to store the population sizes during each generation in the prerun*/
unsigned long popSizeArray[MAXGEN+1];		/*array used to store the population sizes for each gen for a given cycle*/
unsigned long curPopSize;					/*The current population size at any time*/
unsigned long curGenNum;					/*the current generation number*/
unsigned long curRepNum;					/*the current replicate*/
unsigned long sampleSeqs[MAXSAMPLENUM+1];	/*indices of sequences that are to be sampled for the current replicate*/		

/*arrays of flags for processing data*/
char curSeqHitSite[MAXSEQLEN+1]; 			/*1=site already hit in current mutagenesis of a gene, 0=not*/
char curGenHitSite[MAXSEQLEN+1];			/*1=site already hit in current generation (across genes) 0=not*/
char polySite[MAXSEQLEN+1];					/*1=site is segregating in population, 0=not*/		
char crossHitSite[MAXSEQLEN+1];				/*1=site already selected as crossover site in a given generation, 0=not*/

/*variables for data collection from each replicate. */
unsigned long sampleFix01[3+1];				/*Fixations of 0->1 due to sampling (fixed in sample, but not in pop) */
unsigned long sampleFix10[3+1];				/*Fixations of 1->0 due to sampling*/
unsigned long sampleFixNum;					/*Total number of fixations due to sampling */
unsigned long sampleFixPos[MAXSEQLEN+1];	/*Position of fixed (due to sampling) sites in the sample*/
unsigned long sampleSeg01[3+1];				/*segregating sites in the sample with ancSite 0*/
unsigned long sampleSeg10[3+1];				/*segregating sites in the sample with ancSite 1*/
unsigned long sampleSegNum;					/*Total number of segregating sites in the sample*/
unsigned long sampleSegPos[MAXSEQLEN+1];	/*Position of segregating sites in the sample*/

unsigned long popRecNum;					/*Total number of crossovers in the population across all sites */
											
											/*for every fixation in the population, store the derived state, the position in the sequence, and the fitness effect */
unsigned long popFixNum;					/*number of mutations fixed in the replicate - used to access the following arrays*/
char fixSite[MAXFIXATIONS+1];				/*the site that got fixed ie. if 0 was fixed then 0 else 1*/
unsigned long fixSitePos[MAXFIXATIONS+1];	/*Positions of fixed sites in current replicate*/
double fixSiteSval[MAXFIXATIONS+1];			/*S values of sites fixed in current replicate*/

unsigned long repTime;						/*time in seconds for the current replicate*/

unsigned long genMultiHitCount;				/*Total of multiple hits in a generation across all generations for current rep*/
unsigned long segMultiHitCount;				/*Total of multiple hits at segregating sites across all generations*/
/*---------------------------------------------------------------------------------------*/
/*integer random number generator*/
#define getRandLong_(min,max)	((min)+(unsigned long)((1.0-ran2())*((max)-(min)+1)))	/* careful about which ran2 is being used !!!!!!!!!!!!! */
/*---------------------------------------------------------------------------------------*/
/*function declarations*/
int getCtlFile(void);
FILE* fileOpen(char *path,char *mode);
char* getFileString(FILE *fp,char *string,unsigned long maxStringLength);
int readDataWithComments(FILE* fp,char*format,void* target);
int checkParameters(void);
int writeFsimInfo(void);
int initData(void);
void* memAlloc(size_t nmemb,size_t membsize,char* membname);
int doFsim(unsigned long genNum,int flagReset,int repType);
int getPopScenario(int repType,unsigned long genNum);
int resetData(void);
int doMutations(void);
int getHitSites(char origSite,unsigned long baseCount,unsigned long curSeq,unsigned long codPos,unsigned long* hitSites,unsigned long numHits);
unsigned long getPosNthChar(char *site,int codPos,unsigned long hitCod,char origSite);
int doRecombination(void);
int compareUL(const void*a,const void*b);
int getNextGenSeqs(void);
void multdev(double inn[],unsigned long k, unsigned long n, unsigned long nn[]);
double bnldev(double pp,unsigned long n);
int updateCounts(void);
double sampleSval(unsigned long sitePos,char ancSite);
int writeRepData(void);
int getRepData(void);
int getSample(void);
int endSim(void);
int checkBaseCounts(int pos);
double ran2(void);
double poidev(double xm);
double gammln(double xx);
long getTime(int);
/*---------------------------------------------------------------------------------------*/
int main()
{
	/*read the information from control file*/
	getCtlFile();
	/*Initialize data*/
	initData();
	/*writes the simulation information as header of output file
	also outputs the information to the screen*/
	writeFsimInfo();
	/*initrun for system to reach equilibrium,
	DO_RESET - reset data being tracked. Data being tracked include
	counts of fixations, polymorphisms, multihits*/
	printf("Initrun for %lu generation\n",initRunGen);
	doFsim(initRunGen,DO_RESET,INIT_RUN);
	/*initialize timer*/
	getTime(0);
	/*loop for repNum*/
	for(curRepNum=1;curRepNum<=repNum;curRepNum++)
	{
		printf("Replicate %lu of %lu\n",curRepNum,repNum);
		/*prerun for independence between replicates
		reset counts of data being tracked after the run - DO_RESET*/
		printf("PreRun for %lu generations\n",preRunGen); 
		doFsim(preRunGen,DO_RESET,PRE_RUN);
		/*replicate run during which data is collected
		dont reset counts of data being tracked- DONT_RESET
		the data are reset at the end of the prerun, don't reset refers to after the replicate
		reset should not be done prior to writerepdata */ 
		printf("RepRun for %lu generations\n",repRunGen);
		doFsim(repRunGen,DONT_RESET,REP_RUN);
		/*collect data for the current replicate*/
		writeRepData();
		printf("Replicate %lu done in %lus\n",curRepNum,repTime);
	}
	/*free memory and close files*/
	endSim();
	printf("\nfinitho\n");
	getTime(1);
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Reads the input parameters from the control file*/
int getCtlFile(void)
{	
	FILE *fpCtl;
	printf("Enter path+name of control file: ");
	scanf("%s",ctlFileName);
	fpCtl=fileOpen(ctlFileName,"rb");
	/*get input parameters in the following order*/
	getFileString(fpCtl, comments, MAXLINECHAR);
	readDataWithComments(fpCtl,"%s",version);
	readDataWithComments(fpCtl,"%s",fsimVersion);
	readDataWithComments(fpCtl,"%s",machineName);
	readDataWithComments(fpCtl,"%s",outFileName);
	readDataWithComments(fpCtl,"%s",outFileExtn);
	readDataWithComments(fpCtl,"%d",&noMultiHit);
	readDataWithComments(fpCtl,"%lu",&repNum);
	readDataWithComments(fpCtl,"%lu",&initRunGen);
	readDataWithComments(fpCtl,"%lu",&preRunGen);
	readDataWithComments(fpCtl,"%lu",&repRunGen);
	readDataWithComments(fpCtl,"%lu",&initSeqNum);
	readDataWithComments(fpCtl,"%lu",&seqLen);
	readDataWithComments(fpCtl,"%lu",&sampleNum);
	readDataWithComments(fpCtl,"%lf",&c);
	readDataWithComments(fpCtl,"%lf",&u10);
	readDataWithComments(fpCtl,"%lf",&u01);
	/*process filename*/
	strcat(outFileName,"_");
	strcat(outFileName,version);
	strcat(outFileName,".");
	strcat(outFileName,outFileExtn);
	checkParameters();
	fclose(fpCtl);
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Safe file open function. can check for OS specific length issues with filenames
for mac the length is 31 and for OSX it is 255
path - the path to the file+filename or just filename if in same folder
mode - the fopen modes
*/
FILE* fileOpen(char *path,char *mode)
{
    FILE *fp;
	int i,j;
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
		errorOut(("Filename %s exceeded the maximum limit set by the operating system",&path[i]));
    fp=fopen(path,mode);
    if(fp)
        return fp;
    else
        errorOut(("Unable to open %s in %s mode",path,mode));
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
char* getFileString(FILE *fp,char *string,unsigned long maxStringLength)
{
	unsigned long curLength;
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
int readDataWithComments(FILE* fp,char*format,void* target)
{
	char tempc,foundeof=0;
	/*reads the data*/
    if(fscanf(fp,format,target)!=1)
        errorOut(("Corrupted input file or invalid format string(%s)",format));
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
/*Checks the input parameters*/
int checkParameters(void)
{	
	unsigned long maxGen;					/*maximum length of any interval*/
	unsigned long expFixations;				/*expected fixations during the max interval*/
	unsigned long fourZigmaRange;			/*four zigma range for the number of fixations*/
	double maxU;							/*maximum of the two mutation rates*/
	int i,j,imax;							/*temporary ints for looping*/
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
	j=strlen(version);
	imax=strlen(ctlFileName)-j+1;
	i=0;
	while(i<imax)
	{
		if(!strncmp(version,&ctlFileName[i],j))
			break;
		i++;
	}
	if(i==imax)
		errorOut(("version number in control file name doesnt match version number inside file(%s)\n",version));
	/*check if output file already exists*/
	fpOut=fopen(outFileName,"rb");
	if(fpOut)
		errorOut(("output file '%s' exists. delete it or rename it",outFileName));
	else
		fclose(fpOut);
	/*population size*/
	if(initSeqNum>MAXSEQNUM||initSeqNum<MINSEQNUM)
		errorOut(("initSeqNum(%lu) should be between MAXSEQNUM(%lu) and MINSEQNUM(%lu)",initSeqNum,(unsigned long)MAXSEQNUM,(unsigned long)MINSEQNUM));
	if(initSeqNum%2)
		errorOut(("initSeqNum(%lu) has to be a multiple of 2",initSeqNum));
	if(MAXSEQNUM%2)
		errorOut(("MAXSEQNUM(%lu) has to be a multiple of 2",(unsigned long)MAXSEQNUM));
	if(MINSEQNUM%2)
		errorOut(("MAXSEQNUM(%lu) has to be a multiple of 2",(unsigned long)MINSEQNUM));
	/*sample size*/
	if(sampleNum>MINSEQNUM)
		errorOut(("sampleSize(%lu) cannot be greater than MINPOPSIZE(%lu)",sampleNum,(unsigned long)MINSEQNUM));
	if(sampleNum>MAXSAMPLENUM)
		errorOut(("sampleSize(%lu) cannot be greater than MAXSAMPLENUM(%lu)",sampleNum,(unsigned long)MAXSAMPLENUM));
	/*generation sizes*/
	if(initRunGen>MAXGEN||preRunGen>MAXGEN||repRunGen>MAXGEN)
		errorOut(("initRunGen, preRunGen and repRunGen should be smaller than MAXGEN(%lu)",(unsigned long)MAXGEN));
	/*sequence lengths*/
	if(seqLen>MAXSEQLEN)
		errorOut(("seqLen(%lu) cannot be greater than MAXSEQLEN(%lu)",seqLen,(unsigned long)MAXSEQLEN));
	if(seqLen%3)
		errorOut(("seqLen(%lu) has to be a multiple of 3",seqLen));
	if(MAXSEQLEN%3)
		errorOut(("MAXSEQLEN(%lu) has to be a multiple of 3",(unsigned long)MAXSEQLEN));
	/*fixations per interval*/
	maxU=(u10>u01)?u10:u01;
	maxGen=(initRunGen>preRunGen)?initRunGen:preRunGen;
	maxGen=(repRunGen>maxGen)?repRunGen:maxGen;
	expFixations=maxGen*maxU*seqLen;
	fourZigmaRange=expFixations+4*sqrt(expFixations);
	if(MAXFIXATIONS<fourZigmaRange)
		errorOut(("MAXFIXATIONS(%lu) less than four Zigma Range(%lu)",(unsigned long)MAXFIXATIONS,fourZigmaRange));
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*writes the simulation information as the main header in the output file*/
int writeFsimInfo(void)
{
	fprintf(fpOut,"fileName           = %s\n",		outFileName		);
	fprintf(fpOut,"version            = %s\n",		version			);
	fprintf(fpOut,"fsimVersion        = %s\n",		fsimVersion		);
	fprintf(fpOut,"machineName        = %s\n",		machineName		);
	fprintf(fpOut,"MAXSEQNUM          = %lu\n",		MAXSEQNUM		);
	fprintf(fpOut,"MINSEQNUM          = %lu\n",		MINSEQNUM		);
	fprintf(fpOut,"POP_MODEL          = %s\n",		POP_MODEL		);
	fprintf(fpOut,"initSeqNum         = %lu\n",		initSeqNum		);
	fprintf(fpOut,"sampleNum          = %lu\n",		sampleNum		);
	fprintf(fpOut,"seqLen             = %lu\n",		seqLen			);
	fprintf(fpOut,"initRunGen         = %lu\n",		initRunGen		);
	fprintf(fpOut,"preRunGen          = %lu\n",		preRunGen		);
	fprintf(fpOut,"repRunGen          = %lu\n",		repRunGen		);
	fprintf(fpOut,"repNum             = %lu\n",		repNum			);
	fprintf(fpOut,"u01                = %.10f\n",	u01				);
	fprintf(fpOut,"u10                = %.10f\n",	u10				);
	fprintf(fpOut,"c                  = %.10f\n",	c				);
	fprintf(fpOut,"noMultiHit         = %d\n",		noMultiHit		);
	fprintf(fpOut,"S_MODEL            = %s\n",		S_MODEL			);
	fprintf(fpOut,"idum               = %ld\n",		IDUMINIT		);
	fprintf(fpOut,"fsimComments       = %s\n", 		comments		);
	fprintf(fpOut,"//\n");
	fflush(fpOut);
	printf("Writing simulation information to output file\n");
	printf("fileName           = %s\n",				outFileName		);
	printf("version            = %s\n",				version			);
	printf("fsimVersion        = %s\n",				fsimVersion		);
	printf("machineName        = %s\n",				machineName		);
	printf("MAXSEQNUM          = %lu\n",			MAXSEQNUM		);
	printf("MINSEQNUM          = %lu\n",			MINSEQNUM		);
	printf("POP_MODEL          = %s\n",				POP_MODEL		);
	printf("initSeqNum         = %lu\n",			initSeqNum		);
	printf("sampleNum          = %lu\n",			sampleNum		);
	printf("seqLen             = %lu\n",			seqLen			);
	printf("initRunGen         = %lu\n",			initRunGen		);
	printf("preRunGen          = %lu\n",			preRunGen		);
	printf("repRunGen          = %lu\n",			repRunGen		);
	printf("repNum             = %lu\n",			repNum			);
	printf("u01                = %.10f\n",			u01				);
	printf("u10                = %.10f\n",			u10				);
	printf("c                  = %.10f\n",			c				);
	printf("noMultiHit         = %d\n",				noMultiHit		);
	printf("S_MODEL            = %s\n",				S_MODEL			);
	printf("idum               = %ld\n",			IDUMINIT		);
	printf("fsimComments       = %s\n", 			comments		);
	halt;
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*function that initializes all the variables. whenever this function is called
the whole data will be reset. used only at the begining of the whole simulation (not for each replicate) 
some of the variables are reset between prerun and replicate */
int initData(void)
{
	unsigned long numOnes;					/*the number of ones in a sequence*/
	unsigned long numZeros;					/*the number of zeros in a sequence*/
	unsigned long curSite;					/*refers to the current site*/
	unsigned long curSeq;					/*refers to the current sequence*/
	unsigned long randSite;					/*site that is randomly picked to be filled with a 1*/
	double mcu;								/*mcu for the sequence*/
	unsigned long site[MAXSEQLEN+1];		/*temporary array to initialize the sequence*/
	
	printf("Initializing Data\n");
	/*allocate memory for maximum sequences possible*/
	seq=memAlloc(MAXSEQNUM+1,sizeof(struct seq_str),"seq array");
	newSeq=memAlloc(MAXSEQNUM+1,sizeof(struct seq_str),"new seq array");
	
	/*open file*/
	fpOut=fileOpen(outFileName,"wb");
	
	/*finds the initial number of 1s and 0s*/
	mcu=INITMCU;								/* this should not be initialized to 0.50 for MCU model !!!!!!!!!!!!!!!!!!!!!! */
	numOnes=mcu*seqLen;
	numZeros=seqLen-numOnes;
	/*iniitializes the temporary array and fills the first sequence with all zeros
	also initializes the siteData structure*/
	for(curSite=1;curSite<=seqLen;curSite++)
	{
		site[curSite]=curSite;
		seq[1].site[curSite]=0;
	}
	/*find sites with ones*/
	for(curSite=1;curSite<=numOnes;curSite++)
	{
		/*randomly picks a number from those that are not already picked and then
		write a 1 in the seq data structure at the corresponding position*/
		randSite=getRandLong_(curSite,seqLen);
		seq[1].site[site[randSite]]=1;
		site[randSite]=site[curSite];
	}
	/*Initializes the basecounts in the sequence structure and also the sitedata
	structure*/
	seq[1].baseCount0[1]=seq[1].baseCount0[2]=seq[1].baseCount0[3]=0;
	seq[1].baseCount1[1]=seq[1].baseCount1[2]=seq[1].baseCount1[3]=0;
	for(curSite=1;curSite<=seqLen;curSite++)
	{
		/*For each site increment baseCounts, fix siteData elements*/
		if(seq[1].site[curSite]==1)
		{
			seq[1].baseCount1[sitePosToCodPos(curSite)]++;
			siteData[curSite].ancSite=1;
			siteData[curSite].freqs1=initSeqNum;
			siteData[curSite].freqs0=0;
		}
		else
		{
			seq[1].baseCount0[sitePosToCodPos(curSite)]++;
			siteData[curSite].ancSite=0;
			siteData[curSite].freqs0=initSeqNum;
			siteData[curSite].freqs1=0;
		}
		/*Reset Counts in siteData structure*/
		siteData[curSite].fixNum01=0;
		siteData[curSite].fixNum10=0;
		siteData[curSite].mutNum01=0;
		siteData[curSite].mutNum10=0;
		siteData[curSite].recNum=0;
	}
	/*check if everything was ok*/
	if(seq[1].baseCount0[1]+seq[1].baseCount0[2]+seq[1].baseCount0[3]!=numZeros)
		errorOut(("Invalid 0 counts(!=%lu) while initializing sequence",numZeros));
	if(seq[1].baseCount1[1]+seq[1].baseCount1[2]+seq[1].baseCount1[3]!=numOnes)
		errorOut(("Invalid 1 counts(!=%lu) while initializing sequence",numOnes));
	/*copy seq[1] to all sequences*/
	for(curSeq=2;curSeq<=MAXSEQNUM;curSeq++)
	{
		seq[curSeq]=seq[1];
	}
	/*Initializes the sArray*/
	for(curSite=1;curSite<=MAXSEQLEN;curSite++)
	{
		/*sArray[curSite]=sampleSval(curSite,seq[1].site[curSite]);*/
		sArray[curSite]=0.0;
	}
	for(curSite=1;curSite<=MAXGEN+1;curSite++)
	{
		popSizeArray[curSite]=initSeqNum;
	}
	curPopSize=initSeqNum;

	for(curSeq=1;curSeq<=MAXSAMPLENUM;curSeq++)
		sampleSeqs[curSeq]=0;
	
	/*Initializes the flags*/
	for(curSite=1;curSite<=MAXSEQLEN;curSite++)
	{
		curSeqHitSite[curSite]=0;
		curGenHitSite[curSite]=0;
		polySite[curSite]=0;
		crossHitSite[curSite]=0;
	}
	/*Initializes the counts*/
	sampleFix01[1]=sampleFix01[2]=sampleFix01[3]=0;
	sampleFix10[1]=sampleFix10[2]=sampleFix10[3]=0;
	sampleFixNum=0;
	sampleSeg01[1]=sampleSeg01[2]=sampleSeg01[3]=0;
	sampleSeg10[1]=sampleSeg10[2]=sampleSeg10[3]=0;
	sampleSegNum=0;
	for(curSite=1;curSite<=MAXSEQLEN;curSite++)
	{
		sampleFixPos[curSite]=0;
		sampleSegPos[curSite]=0;
	}
	popRecNum=0;
	popFixNum=0;
	for(curSite=1;curSite<=MAXFIXATIONS;curSite++)
	{
		fixSite[curSite]=0;
		fixSitePos[curSite]=0;
		fixSiteSval[curSite]=0;
	}
	repTime=0;
	genMultiHitCount=0;
	segMultiHitCount=0;
	
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*safe memory allocation function. tries to allocate memory and if error then 
print error with variable name and exit
nmemb - the number of elements to be allocated
memmsize - size of each element
membname - name of the variable for which memory is being allocated for errorout
*/
void* memAlloc(size_t nmemb,size_t membsize,char* membname)
{
    void* ptr;
    static unsigned long totmem=0;
    if(!nmemb&&!membsize)
        return &totmem;
    ptr=calloc(nmemb,membsize);
    if(ptr)
    {
        totmem+=nmemb*membsize;
        return ptr;
    }
    else
        errorOut(("Unable to allocate %ld bytes for %s",nmemb*membsize,membname));
    return ptr;
}
/*--------------------------------------------------------------------------------*/
/*the main alogorithm involved in the fsim: mutation, recombination, reproduction
coutfix  are done in that order for the required number of generations
genNum - the number of generations to be run
flagReset - flag which says whether to reset the counts for the current call of fsim
repType - the flag which differentiates between the three types of calls - initrun, prerun and reprun*/
int doFsim(unsigned long genNum,int flagReset,int repType)
{	
	getPopScenario(repType,genNum);
	for(curGenNum=1;curGenNum<=genNum;curGenNum++)
	{
#ifdef DEBUGMODE
		printf("Generation %lu\n",curGenNum);
		halt;
#endif
		curPopSize=popSizeArray[curGenNum];
		doMutations();
#ifdef DEBUGMODE
		checkBaseCounts(1);
#endif
		doRecombination();
#ifdef DEBUGMODE
		checkBaseCounts(2);
#endif
		getNextGenSeqs();
		updateCounts();
#ifdef DEBUGMODE
		checkBaseCounts(3);
#endif
		if(!(curGenNum%500))
			printf("Generation %lu of %lu done\n",curGenNum,genNum);
	}
	if(flagReset==DO_RESET)
		resetData();
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Function that fills the popSizeArray for the current cycle. called from
doFsim at the start of a initrun, prerun or reprun
*/
int getPopScenario(int repType,unsigned long genNum)
{
	unsigned long curPopArrayCell;
#ifdef VERBOSE
	printf("Getting population scenario\n");
#endif
	switch(repType)
	{
		/*InitRun*/
		case 0:	
		{
			for(curPopArrayCell=1;curPopArrayCell<=genNum;curPopArrayCell++)
			{
				popSizeArray[curPopArrayCell]=initSeqNum;
			}
			break;
		}
		/*PreRun*/
		case 1:
		{
			for(curPopArrayCell=1;curPopArrayCell<=genNum;curPopArrayCell++)
			{
				popSizeArray[curPopArrayCell]=initSeqNum;
			}
			/*if it is a prerun then store the population scenario into another array
			so that it can be output after the reprun*/
			for(curPopArrayCell=1;curPopArrayCell<=genNum;curPopArrayCell++)
			{
				preRunPopFluct[curPopArrayCell]=popSizeArray[curPopArrayCell];
			}
			break;
		}
		/*RepRun*/
		case 2:
		{
			for(curPopArrayCell=1;curPopArrayCell<=genNum;curPopArrayCell++)
			{
				popSizeArray[curPopArrayCell]=initSeqNum;
			}
			break;
		}
	}
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Function that resets counts and sitedata after every cycle where data is 
not collected. this function is therefore called after initrun and prerun*/
int resetData(void)
{
	unsigned long curSite;					/*refers to the current site*/
	popFixNum=0;
	popRecNum=0;
	genMultiHitCount=0;
	segMultiHitCount=0;
	for(curSite=1;curSite<=seqLen;curSite++)
	{
		siteData[curSite].fixNum01=0;
		siteData[curSite].fixNum10=0;
		siteData[curSite].mutNum01=0;
		siteData[curSite].mutNum10=0;
		siteData[curSite].recNum=0;
	}
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Function that goes through all the sequences finds the expected number of mutations
and then implements the changes. Uses and updates the polySite array. Resets and uses
the curSeqHitSite array and curGenHitSite array. curGenHitSite array is reset at the 
start of loop and curSeqHitSite array is reset after it is used. So curSeqHitSite array
has to be initialized in the init function. Each mutation results in updating the counts
in the sequence structure . */
int doMutations(void)
{
	int curCodPos;							/*current codon position, 1,2=rep,3=sil*/	
	unsigned long curSite;					/*variable to loop through all flags, also used to implement mutation*/
	unsigned long curSeq;					/*current sequence being mutated*/
	unsigned long curHitSite;				/*Index of array storing hitsites*/
	double expHits10;						/*expected number of hits at 1 sites for a given codpos*/
	double expHits01;						/*expected number of hits at 0 sites for a given codpos*/
	unsigned long numHits10;				/*poisson deviate of the expected number for 1 sites*/
	unsigned long numHits01;				/*poisson deviate of the expected number for 0 sites*/
	unsigned long numChanges10;				/*effective number of changes from 1 to 0*/
	unsigned long hitSites10[MAXSEQLEN+1];	/*indices of 1 sites to be mutated*/
	unsigned long hitSites01[MAXSEQLEN+1];	/*indices of 0 sites to be mutated*/
	
#ifdef DEBUGMODE
	printf("Doing Mutations\n");
	halt;
#endif

	/*reset flags for curGenHitSite. curSeqHitsite are reset in getHitSites function*/
	for(curSite=1;curSite<=seqLen;curSite++)
		curGenHitSite[curSite]=0;
	/*Loop over all sequences*/
	for(curSeq=1;curSeq<=curPopSize;curSeq++)
	{
		/*Loop over positions 1,2,3 of codons*/
		for(curCodPos=1;curCodPos<=3;curCodPos++)
		{	
			/*get expected number of hits for current cod pos*/
			expHits10=u10*seq[curSeq].baseCount1[curCodPos];
			expHits01=u01*seq[curSeq].baseCount0[curCodPos];
			/*take poisson deviates of expected numbers*/
			numHits10=poidev(expHits10);
			numHits01=poidev(expHits01);
			/*limit the poisson deviates to the maximum possible*/
			if(numHits10>seq[curSeq].baseCount1[curCodPos])
				numHits10=seq[curSeq].baseCount1[curCodPos];
			if(numHits01>seq[curSeq].baseCount0[curCodPos])
				numHits01=seq[curSeq].baseCount0[curCodPos];
			/*pick sites at which to change the states*/
			if(numHits10!=0)
			{
				getHitSites(1,seq[curSeq].baseCount1[curCodPos],curSeq,curCodPos,hitSites10,numHits10);
#ifdef VERBOSE
				printf("10 curSeq=%lu curCodPos=%lu\n",curSeq,curCodPos);
				printf("baseCount1 = %lu baseCount0 = %lu\n",seq[curSeq].baseCount1[curCodPos],seq[curSeq].baseCount0[curCodPos]);
#endif
			}
			if(numHits01!=0)
			{
				getHitSites(0,seq[curSeq].baseCount0[curCodPos],curSeq,curCodPos,hitSites01,numHits01);
#ifdef VERBOSE
				printf("01 curSeq=%lu curCodPos=%lu\n",curSeq,curCodPos);
				printf("baseCount1 = %lu baseCount0 = %lu\n",seq[curSeq].baseCount1[curCodPos],seq[curSeq].baseCount0[curCodPos]);
#endif
			}
			/*change sites and update counts for 1->0 mutations also update siteData mutNum10*/
			for(curHitSite=1;curHitSite<=numHits10;curHitSite++)
			{
				curSite=hitSites10[curHitSite];
#ifdef VERBOSE
				printf("%lu\n",curSite);
#endif
#ifdef DEBUGMODE
				if(seq[curSeq].site[curSite]!=1)
					printf("10 not 1\n");
#endif
				seq[curSeq].site[curSite]=0;
				siteData[curSite].mutNum10++;
				siteData[curSite].freqs1--;
				siteData[curSite].freqs0++;
				/* sample new sVal for the sArray for each newly arisen mutation */
				sArray[curSite]=sampleSval(curSite,1);
			}
			/*change sites and update counts for 0->1 mutations also update siteData mutNum01*/
			for(curHitSite=1;curHitSite<=numHits01;curHitSite++)
			{
				curSite=hitSites01[curHitSite];
#ifdef VERBOSE
				printf("%lu\n",curSite);
#endif
#ifdef DEBUGMODE
				if(seq[curSeq].site[curSite]!=0)
					printf("01 not 0\n");
#endif
				seq[curSeq].site[curSite]=1;
				siteData[curSite].mutNum01++;
				siteData[curSite].freqs0--;
				siteData[curSite].freqs1++;
				/*get new sVal for the sArray*/
				sArray[curSite]=sampleSval(curSite,0);
			}
			numChanges10=numHits10-numHits01;
			seq[curSeq].baseCount1[curCodPos]-=numChanges10;
			seq[curSeq].baseCount0[curCodPos]+=numChanges10;
#ifdef VERBOSE
			if(numChanges10)
				printf("baseCount1 = %lu baseCount0 = %lu\n",seq[curSeq].baseCount1[curCodPos],seq[curSeq].baseCount0[curCodPos]);
#endif
		}
	}
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Function which searches for the sites that could be mutated. This is where the
multihit restrictions are implemented. If nomultihit flag is on then no multiple hits
within a generation and at segregating sites are permitted. Regardless of nomultihit,
the same site in a given gene cannot be hit twice in a given generation.
The function also checks every 'limit' passes through the loop if any more mutations are possible. 
This is to prevent cases where the flags prevent all mutations to occur
Parameters passed
origSite 	- The kind of base that gets mutated (1 or 0)
codPos		- The codon position 1,2 repsites 3 silsite
baseCount	- Total number of such bases at all codPos sites
curSeq		- The current sequence where the search is made
hitSites	- The array which is filled with positions where mutations should occur
numHits	 	- The number of hits that should occur*/
int getHitSites(char origSite,unsigned long baseCount,unsigned long curSeq,unsigned long codPos,unsigned long* hitSites,unsigned long numHits)
{
	unsigned long limit;					/*The maximum number of iterations through the loop
											before checking if there are actually any sites left
											that can be mutated. If not then error*/
	unsigned long numRandSearch;			/*the number of searches made for the current function call*/
	unsigned long curNumHits;				/*number of hits located at any time*/
	unsigned long hitCod;					/*a mutation is possible at the hitCod codon with
											origSite at codPos*/
	unsigned long curSite;					/*The actual base position for the hit*/
	unsigned long numUnflagged;				/*The number of unflagged */
	/*limit is the number of random searches made before a check is made on
	the sites to see if the remaining mutations are posssible*/
	limit=numHits*MULTIHITFACTOR;
	numRandSearch=0;
	curNumHits=0;
	/*repeat until all hitsites are found*/
	while(curNumHits<numHits)
	{
		/*keep track of the number of random searches*/
		numRandSearch++;
		/*if the number of searches is greater than or equal to a set limit*/
		if(numRandSearch>=limit)
		{
			/*count the number of unflagged sites*/
			numUnflagged=0;
			if(noMultiHit==1)
			{
				/*if noMultihit is on then check all sites if the polysites flag
				is set. this is because curGenHitSites are also flagged as polySite*/
				for(curSite=codPos;curSite<=seqLen;curSite+=3)
					if(seq[curSeq].site[curSite]==origSite)
						if(polySite[curSite]==0)
							numUnflagged++;
			}
			else
			{
				/*if noMultiHit is off then check all sites if curSeqHitSite flag
				is set.*/
				for(curSite=codPos;curSite<=seqLen;curSite+=3)
					if(seq[curSeq].site[curSite]==origSite)
						if(curSeqHitSite[curSite]==0)
							numUnflagged++;
			}
			if(numUnflagged<(numHits-curNumHits))
			{
				printf("Very High mutation rate\n");
				printf("No unflagged %d sites in the current seq\n",(int)origSite);
				printf("Reduce mutation rate or remove nomultihit flag\n");
				exit(1);
			}
			else
			{
				/*never repeat this process again for the current sequence
				set limit to maximum value*/
				limit=ULONG_MAX;
			}
		}
		/*find the site where a mutation can occur*/
		hitCod=getRandLong_(1,baseCount);
		curSite=getPosNthChar(seq[curSeq].site,codPos,hitCod,origSite);
#ifdef DEBUGMODE
		if(seq[curSeq].site[curSite]!=origSite)
			printf("getPosNthChar returned site which is not correct");
#endif
		/*if multiple hits are not permitted then check multiple hit rules
		and if violated go back and get another site*/
		if(noMultiHit==1)
		{
			/*The order in which these checks are made is significant because
			curGenHitSites are always polySites but the reverse is not true*/
			if(curGenHitSite[curSite]==1)
			{
				genMultiHitCount++;
				continue;
			}
			else if(polySite[curSite]==1)
			{
				segMultiHitCount++;
				continue;
			}
		}
		/*No matter what multiple hit preferences are multiple hits in
		the same sequence at same site is not allowed */
		if(curSeqHitSite[curSite]==0)
		{
			/*update counts for the while loop*/
			curNumHits++;
			/*set flags for polySite, curSeqHit and curGenHit*/
			polySite[curSite]=1;
			curSeqHitSite[curSite]=1;
			curGenHitSite[curSite]=1;
			/*store the curSite in the hitSites array*/
			hitSites[curNumHits]=curSite;
		}
	}
	/*resets the curSeqHitSite array. Since it is reset after it is used first
	the array has to be initialized with all zeros in the initdata function*/
	for(curSite=1;curSite<=numHits;curSite++)
	{
		curSeqHitSite[hitSites[curSite]]=0;
	}
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Funtion which searches for the nth 1 or 0 for a given codon position. for example if
codpos=2,origsite=1 and hitCod=20 then the function will return the site position of
the 20th codon with a 1 at its position 2
site 		- the character array storing the given sequence
codPos		- the position in the codon ie 1,2 or 3
origSite	- the function will search for codons with origSite at codPos
hitCod		- the codon which has been hit
*/
unsigned long getPosNthChar(char *site,int codPos,unsigned long hitCod,char origSite)
{
	unsigned long curCod;					/*refers to current codon*/
	unsigned long curSite;					/*refers to current site*/
	curCod=0;
	curSite=codPos;
	/*loop until curCod is the hitCod codon with origSite at the codPos
	The check for curSite is redundant because we already know that there
	are enough codons to satisfy the condition. So this should just serve
	as an error check. The second condition may be removed */
	while(curSite<=seqLen)
	{
		if(site[curSite]==origSite)
		{
			curCod++;
			if(curCod==hitCod)
				break;
		}
		/*increment current site by 3 to look at codPos in the next codon*/
		curSite+=3;
	}
	/*The error check. This may be removed if needed*/
	if(curCod<hitCod)
	{
		printf("looking for %luth %d\n",hitCod,origSite);
		curCod=0;
		for(curSite=codPos;curSite<=seqLen;curSite+=3)
		{
			printf("%d",site[curSite]);
			if(site[curSite]==origSite)
				curCod++;
		}
		printf("found %lu %ds\n",curCod,origSite);
		errorOut(("Something wrong with baseCounts. cursite>seqlen. Check Code\n"));
	}
	return curSite;
}
/*--------------------------------------------------------------------------------*/
/*Function that goes through all the sequences, randomly pairing them and performs
recombination between the paired sequences. Values are changed only when there is 
actual difference between recombining sites. 
*/
int doRecombination(void)
{
	double expCrossOvers;					/*The expected number of crossovers*/		
	unsigned long maxCrossOvers;			/*The maximum number of crossovers = seqLen-1*/
	unsigned long numCrossOvers;			/*The number of crossovers for the current pair*/
	unsigned long seqIndex[MAXSEQNUM+1];	/*Array to store the indices to pair sequences randomly*/
	unsigned long curSeq;					/*Index of the first of the recombining seq in the seqIndex array*/
	unsigned long curSeqPlusOne;			/*Index of the second recombining seq in the seqIndex array*/
	unsigned long curRecSeq1;				/*Sequence number of the current recombining sequence*/
	unsigned long curRecSeq2;				/*Sequence number of the current recombining sequence*/
	unsigned long nextSeq;					/*Variable used to locate the second parent randombly*/
	unsigned long tempL;					/*temporary long variable*/
	unsigned long crossOverSite[MAXSEQLEN+1];/*array to store the positions of crossovers*/
	unsigned long curCrossSiteNum;			/*The number of cross over sites identified so far*/
	unsigned long curCrossSite;				/*currently identified site for crossover*/
	unsigned long curSeg;					/*The current segment which is being exchanged*/
	unsigned long segStart;					/*Site at which the current segment starts*/
	unsigned long segEnd;					/*Site before which the current segment ends*/
	unsigned long curSite;					/*The current site which is being exchanged*/
	int curCodPos;							/*codPos of the current site*/
	
#ifdef DEBUGMODE
	printf("Doing Recombination\n");
	halt;
#endif
	/*expected number of crossovers is the same for all pairs of sequences */
	expCrossOvers=c*seqLen;
	/*The maximum number of crossovers is one less than the sequence length
	so if there are only 2 bases in the sequence the maximum number of 
	crossovers =1*/
	maxCrossOvers=seqLen-1;
	/*seqIndex array is used to pair up sequences. It is initialiezed with
	numbers 1 to curPopSize*/
	for(curSeq=1;curSeq<=curPopSize;curSeq++)
		seqIndex[curSeq]=curSeq;
	/*The indices of the pairs of sequences that are going to be recombined are
	stored at positions 1 and 2 in the seqIndex array*/	
	curSeq=1;
	curSeqPlusOne=2;
	/*repeat loop until the second of the pair is the last cell in the array*/
	while(curSeqPlusOne<=curPopSize)
	{
		/*find an index to be swapped with the index at curSeqPlusOne position in the 
		seqIndex array. Randomly pick a cell from all the remaining indices*/
		nextSeq=getRandLong_(curSeqPlusOne,curPopSize);
		/*swap the index from the newly located site with that at curSeqPlusOne*/
		tempL=seqIndex[curSeqPlusOne];
		seqIndex[curSeqPlusOne]=seqIndex[nextSeq];
		seqIndex[nextSeq]=tempL;
		/*The indices of the random pair of sequences are at curSeq and curSeqPlusOne*/
		curRecSeq1=seqIndex[curSeq];
		curRecSeq2=seqIndex[curSeqPlusOne];
#ifdef VERBOSE
		printf("%lu-recombing %lu and %lu\n",curSeq,curRecSeq1,curRecSeq2);
		halt;
#endif		
		/*update curSeq and curSeqPlusOne to point to the next pair of cells in the seqIndex array*/
		curSeq+=2;
		curSeqPlusOne+=2;
		/*The number of crossovers is the poisson deviate limited by the maxCrossOvers*/
		numCrossOvers=poidev(expCrossOvers);
		if(numCrossOvers>maxCrossOvers)
			numCrossOvers=maxCrossOvers;
		popRecNum+=numCrossOvers;
		/*If there is going to be any crossOvers then process*/	
		if(numCrossOvers!=0)
		{
			/*crossOverSite corresponds to the site that marks the beginning of
			each segement. so the first segment always starts at 1*/
			crossOverSite[1]=1;
			/*generate numCrossOvers sites using crossHitSite array of flags to 
			identify multiple hits. store the sites in crossOverSite array*/
			curCrossSiteNum=1;
			while(curCrossSiteNum<=numCrossOvers)
			{ 
				curCrossSite=getRandLong_(2,seqLen);
				if(crossHitSite[curCrossSite]==0)
				{
					/*curCrossSiteNum is incremented first before storing the value because
					even when there is only one crossover we take it as 1+1+1 to consider
					the beginning and end of the sequence as two crossover sites. So we will
					have crooOverSite=[1,crossSite,seqLen].numCrossOvers will be incremented
					by 2 after flags are reset
					*/
					curCrossSiteNum++;
					crossOverSite[curCrossSiteNum]=curCrossSite;
					crossHitSite[curCrossSite]=1;
				}
#ifdef VERBOSE
				printf("#");
#endif
			}
			numCrossOvers=numCrossOvers+1;
			/*reset the flag array just at the places they were flagged
			also track the recNum in the siteData structure*/
			curCrossSiteNum=1;
			while(curCrossSiteNum<=numCrossOvers)
			{
				crossHitSite[crossOverSite[curCrossSiteNum]]=0;
				siteData[crossOverSite[curCrossSiteNum]].recNum++;
				curCrossSiteNum++;
			}
			for(curSeg=1;curSeg<=MAXSEQLEN;curSeg++)
			{
				if(crossHitSite[curSeg]==1)
				{
					printf("%lu ",curSeg);
					getchar();getchar();
				}
			}
			/*numCrossOvers is incremented by 2 to account for the beginning and end*/
			numCrossOvers=numCrossOvers+1;
			crossOverSite[numCrossOvers]=seqLen+1;
			/*sort the crossOverSite array. first and last positions are already sorted
			so we are only sorting the others. so sort only if more than 1 site other than
			the beginning and end ie only if numCrossOvers>3*/
			if(numCrossOvers>3)
				qsort(&crossOverSite[2],numCrossOvers-1,sizeof(unsigned long),compareUL);
#ifdef VERBOSE
			printf("Found all crossover sites\n");
			halt;
			for(curSeg=1;curSeg<=numCrossOvers;curSeg++)
				printf("%lu\t",crossOverSite[curSeg]);
			ent;
#endif
			/*Go through alternate segments starting from segment 1. Loop until
			curSeg>numCrossOvers-1*/	
			curSeg=1;
			while(curSeg<numCrossOvers)
			{
				/*A segment starts from one crossOverSite to next crossOverSite-1
				this is why the <segEng check in for loop*/
				segStart=crossOverSite[curSeg];
				segEnd=crossOverSite[curSeg+1];
#ifdef VERBOSE
				printf("Exchanging segment %lu-%lu\n",segStart,segEnd-1);
#endif
				for(curSite=segStart;curSite<segEnd;curSite++)
				{
					/*Even for segments that get exchanged, swap sites only if
					they are different*/
					if(polySite[curSite])
					{
						if(seq[curRecSeq1].site[curSite]==1)
						{
							if (seq[curRecSeq2].site[curSite]==0)
							{
								/*site in first seq is 1 and in second is 0
								so swap and update counts*/
								curCodPos=sitePosToCodPos(curSite);
								seq[curRecSeq1].site[curSite]=0;
								seq[curRecSeq1].baseCount0[curCodPos]++;
								seq[curRecSeq1].baseCount1[curCodPos]--;
								seq[curRecSeq2].site[curSite]=1;
								seq[curRecSeq2].baseCount0[curCodPos]--;
								seq[curRecSeq2].baseCount1[curCodPos]++;
							}
						}
						else 
						{
							if (seq[curRecSeq2].site[curSite]==1)
							{
								/*site in first seq is 0 and in second is 1
								so swap and update counts*/
								curCodPos=sitePosToCodPos(curSite);
								seq[curRecSeq1].site[curSite]=1;
								seq[curRecSeq1].baseCount0[curCodPos]--;
								seq[curRecSeq1].baseCount1[curCodPos]++;
								seq[curRecSeq2].site[curSite]=0;
								seq[curRecSeq2].baseCount0[curCodPos]++;
								seq[curRecSeq2].baseCount1[curCodPos]--;
							}
						}
					}
				}
				/*skip one segment and go to the one after that*/
				curSeg=curSeg+2;
			}
		}
	}
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Function to compare to unsigned long values given that they are passed
as constant void pointers. Just subtracts and returns values. To be used
with qsort*/
int compareUL(const void*a,const void*b)
{
	if((*(unsigned long*)a)>(*(unsigned long*)b))
		return 1;
	else if((*(unsigned long*)a)<(*(unsigned long*)b))
		return -1;
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Function that evaluates the fitness of all individuals and then using multinomial 
sampling finds the numbers of offspring for each individual and creates the next
generation of individuals. Function uses sArray to find fitness of individuals*/
int getNextGenSeqs(void)
{
	unsigned long curSeq;					/*refers to current sequence*/
	unsigned long curParent;				/*refers to the current parent sequence*/
	unsigned long curChild;					/*refers to all child sequences of a given parent*/
	unsigned long curSite;					/*refers to current site*/
	unsigned long newPopSize;				/*the new population size*/
	double wData[MAXSEQNUM+1];				/*array that stores the weight for each sequence*/
	double wSum;							/*sum of wData*/
	double freqs[MAXSEQNUM+1];				/*expected frequencies of current individuals*/
	unsigned long newFreqs[MAXSEQNUM+1];	/*obtained frequencies of offsprings for current generation*/
	void *tempPV;							/*temporary pointer variable for swapping sequence pointers*/
	
#ifdef DEBUGMODE
	printf("Finding next gen seqs\n");
	halt;
#endif	
	/*Initialize weights to 1.0 each*/
	for(curSeq=1;curSeq<=curPopSize;curSeq++)
	{
		wData[curSeq]=1.0;
	}
	/*Loops through all variable sites and finds the weights for each individual*/
	for(curSite=1;curSite<=seqLen;curSite++)
	{	
		if(polySite[curSite])
		{
			for(curSeq=1;curSeq<=curPopSize;curSeq++)
			{
				if(seq[curSeq].site[curSite]!=siteData[curSite].ancSite)
				{
					wData[curSeq]*=(1.0+sArray[curSite]);						/* multiplicative fitness implemented !!! */
				}
			}
		}
	}
	/*Finds the sum of weights to get relative weights*/
	wSum=0.0;
	for(curSeq=1;curSeq<=curPopSize;curSeq++)
	{	
		wSum+=wData[curSeq];
	}
	/*expected frequencies in the nextgeneration = relative weights*/
	for(curSeq=1;curSeq<=curPopSize;curSeq++)
	{	
		freqs[curSeq]=wData[curSeq]/wSum;
///        freqs[curSeq]=1.0/curPopSize;
///        printf("%f -> %f\n", wData[curSeq]/wSum, freqs[curSeq]);
	}
	/*the population size for a given generation is given in the popSizeArray*/
	newPopSize=popSizeArray[curGenNum+1];
	/*take multinomial deviates on the expected frequencies*/
	multdev(freqs,curPopSize,newPopSize,newFreqs);
	/*create new generation of individuals from the new frequences from 
	multinomial sampling*/
	curSeq=1;
	for(curParent=1;curParent<=curPopSize;curParent++)
	{
		for(curChild=1;curChild<=newFreqs[curParent];curChild++)
		{
			newSeq[curSeq]=seq[curParent];
			curSeq++;
		}
	}
	/*check if sampling is working correctly*/
	if((curSeq-1)!=newPopSize)
		errorOut(("Invalid number of new sequences(%lu!=%lu) in getNextGenSeqs",curSeq,newPopSize));
	/*population size=population size in the new generation*/	
	curPopSize=newPopSize;
	/*swap pointers*/
	tempPV=seq;
	seq=newSeq;
	newSeq=tempPV;
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Multinomial deviate function. Function fills in the array nn with the number of 
successes for the corresponding probability of successes given in the inn array and
the total number of experiments n. k is the size of the two arrays.(1 based arrays)
Function Created:  Hiroshi Akashi
*/
void multdev(double inn[],unsigned long k, unsigned long n, unsigned long nn[])
{
	unsigned long i,r;
	double sum;
	sum = 0.0;
	for(i=1;i<=k;i++)
	{
		if( (n > 0) && (sum < 1.0) )
		{
			r = (unsigned long)(bnldev(inn[i]/(1.0 - sum), n));
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
double bnldev(double pp,unsigned long n)
{
	unsigned long j;
	static long nold=(-1);
	double am,em,g,angle,p,bnl,sq,t,y;
	static double pold=(-1.0),pc,plog,pclog,en,oldg;
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
                    em=(unsigned long)(em);
                    t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
                            -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
            } while (ran2() > t);
            bnl=em;
    }
    if (p != pp) bnl=n-bnl;
    return bnl;
}
/*--------------------------------------------------------------------------------*/
/*During reproduction sequences may get lost due to multinomial sampling and/or selection.
Polymorphic sites may become monomorphic or the derived state might get fixed
in the generation. This function updates the siteData structure and the polySites
flag. Also stores the fixation positions,svals and ancsite so that these can be 
stored in the outputfile
*/
int updateCounts(void)
{
	unsigned long curSite;					/*refers to the current site*/
	unsigned long curSeq;					/*refers to the current sequence*/
	unsigned long baseCount0;				/*number of 0 bases at curSite across all seqs*/
	unsigned long baseCount1;				/*number of 1 bases at curSite across all seqs*/
	
#ifdef DEBUGMODE
	printf("Updating Counts\n");
	halt;
#endif
	/*Loop across all the sites*/
	for(curSite=1;curSite<=seqLen;curSite++)
	{
		/*check only if the polySite flag is set*/
		if(polySite[curSite]==1)
		{
			/*Find the frequencies of 0s and 1s at the site*/
			baseCount0=0;
			baseCount1=0;
			for(curSeq=1;curSeq<=curPopSize;curSeq++)
			{
				if(seq[curSeq].site[curSite]==0)
					baseCount0++;
				else
					baseCount1++;
			}
			/*update the counts in the siteData structure*/
			siteData[curSite].freqs0=baseCount0;
			siteData[curSite].freqs1=baseCount1;
			/*Check if any of the basecounts is 0*/
			if(baseCount1==0)
			{
				/*if so check if the ancestral site is same as the base whose count is zero
				if so a fixation. here check if 0 has been fixed*/
				if(siteData[curSite].ancSite==1)
				{
					/*if fixation update counts in the siteData structure*/
					siteData[curSite].ancSite=0;
					siteData[curSite].fixNum10++;
					/*store the svalues positions and fixed site to be output in the writerepdata*/
					popFixNum++;
					fixSite[popFixNum]=0;
					fixSitePos[popFixNum]=curSite;
					fixSiteSval[popFixNum]=sArray[curSite];
					/*at every fixation resample for the sArray*/
					sArray[curSite]=0.0;
				}
				polySite[curSite]=0;
			}
			else if(baseCount0==0)
			{
				/*if so check if the ancestral site is same as the base whose count is zero
				if so a fixation. here check if 1 has been fixed*/
				if(siteData[curSite].ancSite==0)
				{
					/*if fixation update counts in the siteData structure*/
					siteData[curSite].ancSite=1;
					siteData[curSite].fixNum01++;
					/*store the svalues positions and fixed site to be output in the writerepdata*/
					popFixNum++;
					fixSite[popFixNum]=1;
					fixSitePos[popFixNum]=curSite;
					fixSiteSval[popFixNum]=sArray[curSite];
					/*at every fixation resample for the sArray*/
					sArray[curSite]=0.0;
				}
				polySite[curSite]=0;
			}
		}
	}
	/*If number of fixations greater than maximum permitted then error*/
	if(popFixNum>MAXFIXATIONS)
		errorOut(("Exceeded maximum fixations permitted, redefine MAXFIXATIONS"));
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Function used to fill and modify the sArray. This function is called
when initializing and when fixations occur*/
double sampleSval(unsigned long sitePos,char ancSite)
{
	int codPos;
	const double Nes = 0.00001;
	const double Ne = 500.0;
	const double s01 = Nes/Ne;
	const double s02 = Nes/Ne;
	const double s03 = Nes/Ne;
	const double s11 = Nes/Ne;
	const double s12 = Nes/Ne;
	const double s13 = Nes/Ne;
	const double prob = 1.0/10.0;					/*the fraction of times the function returns a nonzero value*/
													/*changed the variables to consts - AJ 07/21/03*/
	if(ran2() > prob)
		return 0.0;
	codPos=sitePosToCodPos(sitePos);
	switch(codPos)
	{
		case	1:	
					if(ancSite==1)
						return s11;
					else
						return s01;
					break;
		case	2:	
					if(ancSite==1)
						return s12;
					else
						return s02;
					break;
		case	3:	
					if(ancSite==1)
						return s13;
					else
						return s03;
					break;
	}
	return 0.0;
}
/*--------------------------------------------------------------------------------*/
/*Writes the information from each replicate into the output file. Each repdata
consists of a set of counts and the data from the samples. population fluctuation
information, fixation information including position,fixed site and Svals, segregating
sites information and the common ancestor of the samples are also printed ou*/
int writeRepData(void)
{
	unsigned long i,j;						/*for loop*/
	unsigned long numPopPeriods;			/*number of periods for the population fluctuation scenario*/
	unsigned long startPeriod;				/*starting generation of a population period*/
	unsigned long endPeriod;				/*ending generation of a population period*/
	unsigned long curSeq;					/*refers to the current sampled sequence*/
	
	/*sample sequences and collect data from the samples*/
	getRepData();
	printf("Writing Replicate information to file\n");
	/*Outputs the counts obtained during the current replicate*/
	fprintf(fpOut,"Replicate:%6lu\n",curRepNum);
	fprintf(fpOut,"----------------\n");
	fprintf(fpOut,"sampleFix01[1]     = %lu\n",		sampleFix01[1]	);
	fprintf(fpOut,"sampleFix01[2]     = %lu\n",		sampleFix01[2]	);
	fprintf(fpOut,"sampleFix01[3]     = %lu\n",		sampleFix01[3]	);
	fprintf(fpOut,"sampleFix10[1]     = %lu\n",		sampleFix10[1]	);
	fprintf(fpOut,"sampleFix10[2]     = %lu\n",		sampleFix10[2]	);
	fprintf(fpOut,"sampleFix10[3]     = %lu\n",		sampleFix10[3]	);
	fprintf(fpOut,"sampleSeg01[1]     = %lu\n",		sampleSeg01[1]	);
	fprintf(fpOut,"sampleSeg01[2]     = %lu\n",		sampleSeg01[2]	);
	fprintf(fpOut,"sampleSeg01[3]     = %lu\n",		sampleSeg01[3]	);
	fprintf(fpOut,"sampleSeg10[1]     = %lu\n",		sampleSeg10[1]	);
	fprintf(fpOut,"sampleSeg10[2]     = %lu\n",		sampleSeg10[2]	);
	fprintf(fpOut,"sampleSeg10[3]     = %lu\n",		sampleSeg10[3]	);
	fprintf(fpOut,"popRecNum          = %lu\n",		popRecNum		);
	fprintf(fpOut,"genMultiHitCount   = %lu\n",		genMultiHitCount);
	fprintf(fpOut,"segMultiHitCount   = %lu\n",		segMultiHitCount);
	fprintf(fpOut,"repTime            = %lu\n",		repTime			);
	fprintf(fpOut,"\n");
	/*outputs the population fluctuation scenario for the prerun*/
	fprintf(fpOut,"Population fluctuation for preRun(startGen:endGen:popSize)\n");
	fprintf(fpOut,"----------------------------------------------------------\n");
	/*gets the number of periods*/
	numPopPeriods=1;
	i=2;
	while(i<=preRunGen)
	{
		if(preRunPopFluct[i]!=preRunPopFluct[i-1])
			numPopPeriods++;
		i++;
	}
	fprintf(fpOut,"%lu\n",numPopPeriods);
	/*prints out the population fluctuation for the last period*/
	i=2;
	startPeriod=1;
	numPopPeriods=1;
	while(i<=preRunGen)
	{
		if(preRunPopFluct[i]!=preRunPopFluct[i-1])
		{
			endPeriod=i-1;
			fprintf(fpOut,"%lu:%lu:%lu ",startPeriod,endPeriod,preRunPopFluct[startPeriod]);
			startPeriod=i;
			numPopPeriods++;
		}
		i++;
	}
	if(numPopPeriods==1)
		fprintf(fpOut,"%d:%lu:%lu ",1,preRunGen,preRunPopFluct[1]);
	fprintf(fpOut,"\n");
	fprintf(fpOut,"\n");
	
	/*outputs the population fluctuation scenario*/
	fprintf(fpOut,"Population fluctuation for repRun(startGen:endGen:popSize)\n");
	fprintf(fpOut,"----------------------------------------------------------\n");
	/*gets the number of periods*/
	numPopPeriods=1;
	i=2;
	while(i<=repRunGen)
	{
		if(popSizeArray[i]!=popSizeArray[i-1])
			numPopPeriods++;
		i++;
	}
	fprintf(fpOut,"%lu\n",numPopPeriods);
	/*prints out the population fluctuation for the last period*/
	i=2;
	startPeriod=1;
	numPopPeriods=1;
	while(i<=repRunGen)
	{
		if(popSizeArray[i]!=popSizeArray[i-1])
		{
			endPeriod=i-1;
			fprintf(fpOut,"%lu:%lu:%lu ",startPeriod,endPeriod,popSizeArray[startPeriod]);
			startPeriod=i;
			numPopPeriods++;
		}
		i++;
	}
	if(numPopPeriods==1)
		fprintf(fpOut,"%d:%lu:%lu ",1,repRunGen,popSizeArray[1]);
	fprintf(fpOut,"\n");
	fprintf(fpOut,"\n");
	fprintf(fpOut,"Fixed Svals(fixedSite:fixedSval:pos)\n");
	fprintf(fpOut,"-----------------------------------\n");
	fprintf(fpOut,"%lu\n",popFixNum);
	for(i=1;i<=popFixNum;i++)
	{
		fprintf(fpOut,"%d:%.10f:%lu ",(int)fixSite[i],fixSiteSval[i],fixSitePos[i]);
	}
	fprintf(fpOut,"\n");
	fprintf(fpOut,"\n");
	/*prints out the sitedata information plus sArray*/
	fprintf(fpOut,"SiteData(ancSite:sVal:freqs0:freqs1:fix01:fix10:mut01:mut10:rec)\n");
	fprintf(fpOut,"----------------------------------------------------------------\n");
	for(i=1;i<=seqLen;i++)
	{
		fprintf(fpOut,"%d:",(int)siteData[i].ancSite);
		fprintf(fpOut,"%.10f:",sArray[i]);
		fprintf(fpOut,"%lu:",siteData[i].freqs0);
		fprintf(fpOut,"%lu:",siteData[i].freqs1);
		fprintf(fpOut,"%lu:",siteData[i].fixNum01);
		fprintf(fpOut,"%lu:",siteData[i].fixNum10);
		fprintf(fpOut,"%lu:",siteData[i].mutNum01);
		fprintf(fpOut,"%lu:",siteData[i].mutNum10);
		fprintf(fpOut,"%lu ",siteData[i].recNum);
	}
	fprintf(fpOut,"\n");
	fprintf(fpOut,"\n");
	fprintf(fpOut,"sampleFixPos\n");
	fprintf(fpOut,"------------\n");
	fprintf(fpOut,"%lu\n",sampleFixNum);
	for(j=1;j<=sampleFixNum;j++)
	{
		fprintf(fpOut,"%lu ",sampleFixPos[j]);
	}
	fprintf(fpOut,"\n");
	fprintf(fpOut,"\n");
	fprintf(fpOut,"sampleSegPos\n");
	fprintf(fpOut,"------------\n");
	fprintf(fpOut,"%lu\n",sampleSegNum);
	for(j=1;j<=sampleSegNum;j++)
	{
		fprintf(fpOut,"%lu ",sampleSegPos[j]);
	}
	fprintf(fpOut,"\n");
	fprintf(fpOut,"\n");
	/*prints out the sample sequences*/
	for(i=1;i<=sampleNum;i++)
	{
		curSeq=sampleSeqs[i];
		for(j=1;j<=sampleSegNum;j++)
		{
			fprintf(fpOut,"%d",seq[curSeq].site[sampleSegPos[j]]);
		}
		fprintf(fpOut,"\n");
	}
	fprintf(fpOut,"//\n");
	fflush(fpOut);
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Collects the information from each replicate. This involves picking a sample set
from the population and collecting information from the sample.*/
int getRepData(void)
{
	unsigned long curSite;					/*refers to the current site*/
	unsigned long freq0;					/*frequencies of 0 at the current site*/
	unsigned long freq1;					/*frequencies of 1 at the current site*/
	unsigned long curSample;				/*refers to the current sample num*/

	printf("Getting Replicate Information\n");
	/*reset counts*/
	sampleFix01[1]=sampleFix01[2]=sampleFix01[3]=0;
	sampleFix10[1]=sampleFix10[2]=sampleFix10[3]=0;
	sampleFixNum=0;
	sampleSeg01[1]=sampleSeg01[2]=sampleSeg01[3]=0;
	sampleSeg10[1]=sampleSeg10[2]=sampleSeg10[3]=0;
	sampleSegNum=0;
	/*picks sampleNum sequences*/
	getSample();
	for(curSite=1;curSite<=seqLen;curSite++)
	{
		/*if the site is not polymorphic in the population then no need to process it*/
		if(polySite[curSite]==1)
		{
			freq0=freq1=0;
			for(curSample=1;curSample<=sampleNum;curSample++)
			{
				if(seq[sampleSeqs[curSample]].site[curSite]==0)
					freq0++;
				else
					freq1++;
			}
			/*if site is polymorphic in sample*/
			if((freq1!=0)&&(freq0!=0))
			{
				/*depending on the ancestral state update the counters*/
				if(siteData[curSite].ancSite==1)
				{	
					sampleSeg10[sitePosToCodPos(curSite)]++;
				}
				else
				{
					sampleSeg01[sitePosToCodPos(curSite)]++;
				}
				/*keep track of the segregating sites*/
				sampleSegNum++;
				sampleSegPos[sampleSegNum]=curSite;
			}
			else
			{	
				/*fixation of 0 due to sampling*/
				if((freq1==0)&&(siteData[curSite].ancSite==1))
				{
					sampleFix10[sitePosToCodPos(curSite)]++;
					/*keep track of fixations*/
					sampleFixNum++;
					sampleFixPos[sampleFixNum]=curSite;
				}
				/*fixation of 1 due to sampling*/
				if((freq0==0)&&(siteData[curSite].ancSite==0))
				{
					sampleFix01[sitePosToCodPos(curSite)]++;
					/*keep track of fixations*/
					sampleFixNum++;
					sampleFixPos[sampleFixNum]=curSite;
				}
			}
		}
	}
	repTime=getTime(4);
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*randomly picks the current samplenum sequences from the population*/
int getSample(void)
{
	unsigned long curSample;					/*refers to the current sample*/
	unsigned long randSeq;						/*randomly selected sequence*/
	unsigned long curSeq;						/*loop through all samples*/
	/*repeat loop until all random sequences are found*/
	curSample=0;
	while(curSample<sampleNum)
	{
		/*pick a random sequence*/
		randSeq=getRandLong_(1,curPopSize);
		/*check if the sequence has already been picked*/
		for(curSeq=1;curSeq<=curSample;curSeq++)
		{
			if(sampleSeqs[curSeq]==randSeq)
				break;
		}
		/*if not already picked then add it to the list*/
		if(curSeq>curSample)
		{
			curSample++;
			sampleSeqs[curSample]=randSeq;
		}
	}
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*function for freeing the memory allocated and closing all files*/
int endSim(void)
{
	free(seq);
	free(newSeq);
	fclose(fpOut);
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*function to check basecounts*/
int checkBaseCounts(int pos)
{
	unsigned long i;
	unsigned long curSeq;
	unsigned long curSite;
	unsigned long baseCount0[3+1];
	unsigned long baseCount1[3+1];
	for(curSeq=1;curSeq<=curPopSize;curSeq++)
	{
		baseCount0[1]=baseCount0[2]=baseCount0[3]=0;
		baseCount1[1]=baseCount1[2]=baseCount1[3]=0;
		for(curSite=1;curSite<=seqLen;curSite++)
		{
			if(seq[curSeq].site[curSite]==1)
				baseCount1[(curSite+2)%3+1]++;
			else
				baseCount0[(curSite+2)%3+1]++;
		}
		for(i=1;i<=3;i++)
		{
			if(baseCount0[i]!=seq[curSeq].baseCount0[i])
			{
				printf("Error%d in baseCount0 expected %lu found %lu for codPos %lu in seq %lu\n",pos,baseCount0[i],seq[curSeq].baseCount0[i],i,curSeq);
				halt;
			}	
			if(baseCount1[i]!=seq[curSeq].baseCount1[i])
			{
				printf("Error%d in baseCount1 expected %lu found %lu for codPos %lu in seq %lu\n",pos,baseCount1[i],seq[curSeq].baseCount1[i],i,curSeq);
				halt;
			}	
		}
	}
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*Modified Ran2 function. Original function from the Numerical Recipes in C IIed. 
This function uses a global idum. returns random number in the range (0,1].
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
of the function to get (0,1) output.
*/
/*
#define EPS DBL_EPSILON
#define RNMX (1.0-EPS)
*/

double ran2(void)
{
	int j;
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
	/*the output will be (0,1]*/
	return  (double)AM2*iy;
/*	
	Here the output will be (0,1)
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
       		em=(unsigned long)(em);
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
	int j;

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
long getTime(int m)
{
	static unsigned long start=0,end=0;
	static unsigned long sttemp=0;
	unsigned long diff;
	if(m==0)
	{
		start=time(NULL);
	}
	/*time from start*/
	else if(m==1)
	{
		end=time(NULL);
		diff=difftime(end,start);
		printf("\nTotal Time from start:  %02lu:%02lu:%02lu = %02lu:%02lu\n",diff/3600,(diff%3600)/60,diff%60,diff/60,diff%60);
	}
	/*time since last call*/
	else if(m==2)
	{
		if(end)sttemp=end;
		else sttemp=start;
		end=time(NULL);
		diff=difftime(end,sttemp);
		printf("Last Interval:  %02lu:%02lu:%02lu = %02lu:%02lu\n",diff/3600,(diff%3600)/60,diff%60,diff/60,diff%60);
	}
	/*time from start*/
	else if(m==3)
	{
		end=time(NULL);
		diff=difftime(end,start);
		return diff;
	}
	/*time since last call*/
	else if(m==4)
	{
		if(end)sttemp=end;
		else sttemp=start;
		end=time(NULL);
		diff=difftime(end,sttemp);
		return diff;
	}
	return 0;
}
/*---------------------------------------------------------------------------------------*/
