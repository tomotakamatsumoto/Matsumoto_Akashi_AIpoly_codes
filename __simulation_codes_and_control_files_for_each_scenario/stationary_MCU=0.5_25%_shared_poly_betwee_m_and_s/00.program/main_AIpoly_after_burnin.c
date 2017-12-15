/*
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>				
#include <errno.h>

#include "header1.h"
#include "Userpref_after_burnin.h"
#include "header2_after_burnin.h"
//#include "MT.h"

int check=0;
long population_size;
/*---------------------------------------------------------------------------------------*/
int main(void)
{
	struct userPref_str *userPrefp;
	
	printf("Start of Program fsim_v3 ...\n");
	readControlFile(&userPrefp);	/* Reads the control file and stores the data in userPrefp structure*/

	programTimer(0);				/*start timer*/

//	getCtlFile();

	start_time = clock();
	/*Initialize data*/
	initData(userPrefp);
	
	/* get static fitness array - depends on number of 1's */
	//get_static_w_dat(userPrefp);
	
	/*writes the simulation information as header of output file
	also outputs the information to the screen*/
	writeFsimInfo(userPrefp);
    int tomo;
    printf ("MRCA of initial population is \n");
    for (tomo=0;tomo<=userPrefp->seqLen;tomo++) {
        printf("%c", siteData[tomo].ancSite);
    }
    printf ("\n");

	
	/*loop for repNum*/
	for(curRepNum=1; curRepNum<=userPrefp->repNum; curRepNum++)
	{
        check=1;

		printf ("main run started\n");
		doFsim(userPrefp, userPrefp->repRunGen, DONT_RESET, REP_RUN);
		
		/*collect data for the current replicate*/
		writeRepData(userPrefp, curGenNum);
		//		printf("\tReplicate %ld done in %lus\n", curRepNum, repTime);
	}
    printf ("main run finished\n");
	
	/*free memory and close files*/
	endSim(userPrefp);
	printf("\nfinitho\n");
//	getTime(1);
	end_time = clock();
	getTime_print(userPrefp, 1);
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*
void get_static_w_dat(struct userPref_str *userPrefp)

{
	double s_val = 0;
    double ParaX;
	long i, debug = 1;
	
    cod_w_list_A[0] = userPrefp->f000_a; cod_w_list_A[1] = userPrefp->f001_a; cod_w_list_A[2] = userPrefp->f002_a; cod_w_list_A[3] = userPrefp->f003_a;
    cod_w_list_A[4] = userPrefp->f010_a; cod_w_list_A[5] = userPrefp->f011_a; cod_w_list_A[6] = userPrefp->f012_a; cod_w_list_A[7] = userPrefp->f013_a;
    cod_w_list_A[8] = userPrefp->f020_a; cod_w_list_A[9] = userPrefp->f021_a; cod_w_list_A[10] = userPrefp->f022_a; cod_w_list_A[11] = userPrefp->f023_a;
    cod_w_list_A[12] = userPrefp->f030_a; cod_w_list_A[13] = userPrefp->f031_a; cod_w_list_A[14] = userPrefp->f032_a; cod_w_list_A[15] = userPrefp->f033_a;
    cod_w_list_A[16] = userPrefp->f100_a; cod_w_list_A[17] = userPrefp->f101_a; cod_w_list_A[18] = userPrefp->f102_a; cod_w_list_A[19] = userPrefp->f103_a;
    cod_w_list_A[20] = userPrefp->f110_a; cod_w_list_A[21] = userPrefp->f111_a; cod_w_list_A[22] = userPrefp->f112_a; cod_w_list_A[23] = userPrefp->f113_a;
    cod_w_list_A[24] = userPrefp->f120_a; cod_w_list_A[25] = userPrefp->f121_a; cod_w_list_A[26] = userPrefp->f122_a; cod_w_list_A[27] = userPrefp->f123_a;
    cod_w_list_A[28] = userPrefp->f130_a; cod_w_list_A[29] = userPrefp->f131_a; cod_w_list_A[30] = userPrefp->f132_a; cod_w_list_A[31] = userPrefp->f133_a;
    cod_w_list_A[32] = userPrefp->f200_a; cod_w_list_A[33] = userPrefp->f201_a; cod_w_list_A[34] = userPrefp->f202_a; cod_w_list_A[35] = userPrefp->f203_a;
    cod_w_list_A[36] = userPrefp->f210_a; cod_w_list_A[37] = userPrefp->f211_a; cod_w_list_A[38] = userPrefp->f212_a; cod_w_list_A[39] = userPrefp->f213_a;
    cod_w_list_A[40] = userPrefp->f220_a; cod_w_list_A[41] = userPrefp->f221_a; cod_w_list_A[42] = userPrefp->f222_a; cod_w_list_A[43] = userPrefp->f223_a;
    cod_w_list_A[44] = userPrefp->f230_a; cod_w_list_A[45] = userPrefp->f231_a; cod_w_list_A[46] = userPrefp->f232_a; cod_w_list_A[47] = userPrefp->f233_a;
    cod_w_list_A[48] = userPrefp->f300_a; cod_w_list_A[49] = userPrefp->f301_a; cod_w_list_A[50] = userPrefp->f302_a; cod_w_list_A[51] = userPrefp->f303_a;
    cod_w_list_A[52] = userPrefp->f310_a; cod_w_list_A[53] = userPrefp->f311_a; cod_w_list_A[54] = userPrefp->f312_a; cod_w_list_A[55] = userPrefp->f313_a;
    cod_w_list_A[56] = userPrefp->f320_a; cod_w_list_A[57] = userPrefp->f321_a; cod_w_list_A[58] = userPrefp->f322_a; cod_w_list_A[59] = userPrefp->f323_a;
    cod_w_list_A[60] = userPrefp->f330_a; cod_w_list_A[61] = userPrefp->f331_a; cod_w_list_A[62] = userPrefp->f332_a; cod_w_list_A[63] = userPrefp->f333_a;
    
    cod_w_list_B[0] = userPrefp->f000_b; cod_w_list_B[1] = userPrefp->f001_b; cod_w_list_B[2] = userPrefp->f002_b; cod_w_list_B[3] = userPrefp->f003_b;
    cod_w_list_B[4] = userPrefp->f010_b; cod_w_list_B[5] = userPrefp->f011_b; cod_w_list_B[6] = userPrefp->f012_b; cod_w_list_B[7] = userPrefp->f013_b;
    cod_w_list_B[8] = userPrefp->f020_b; cod_w_list_B[9] = userPrefp->f021_b; cod_w_list_B[10] = userPrefp->f022_b; cod_w_list_B[11] = userPrefp->f023_b;
    cod_w_list_B[12] = userPrefp->f030_b; cod_w_list_B[13] = userPrefp->f031_b; cod_w_list_B[14] = userPrefp->f032_b; cod_w_list_B[15] = userPrefp->f033_b;
    cod_w_list_B[16] = userPrefp->f100_b; cod_w_list_B[17] = userPrefp->f101_b; cod_w_list_B[18] = userPrefp->f102_b; cod_w_list_B[19] = userPrefp->f103_b;
    cod_w_list_B[20] = userPrefp->f110_b; cod_w_list_B[21] = userPrefp->f111_b; cod_w_list_B[22] = userPrefp->f112_b; cod_w_list_B[23] = userPrefp->f113_b;
    cod_w_list_B[24] = userPrefp->f120_b; cod_w_list_B[25] = userPrefp->f121_b; cod_w_list_B[26] = userPrefp->f122_b; cod_w_list_B[27] = userPrefp->f123_b;
    cod_w_list_B[28] = userPrefp->f130_b; cod_w_list_B[29] = userPrefp->f131_b; cod_w_list_B[30] = userPrefp->f132_b; cod_w_list_B[31] = userPrefp->f133_b;
    cod_w_list_B[32] = userPrefp->f200_b; cod_w_list_B[33] = userPrefp->f201_b; cod_w_list_B[34] = userPrefp->f202_b; cod_w_list_B[35] = userPrefp->f203_b;
    cod_w_list_B[36] = userPrefp->f210_b; cod_w_list_B[37] = userPrefp->f211_b; cod_w_list_B[38] = userPrefp->f212_b; cod_w_list_B[39] = userPrefp->f213_b;
    cod_w_list_B[40] = userPrefp->f220_b; cod_w_list_B[41] = userPrefp->f221_b; cod_w_list_B[42] = userPrefp->f222_b; cod_w_list_B[43] = userPrefp->f223_b;
    cod_w_list_B[44] = userPrefp->f230_b; cod_w_list_B[45] = userPrefp->f231_b; cod_w_list_B[46] = userPrefp->f232_b; cod_w_list_B[47] = userPrefp->f233_b;
    cod_w_list_B[48] = userPrefp->f300_b; cod_w_list_B[49] = userPrefp->f301_b; cod_w_list_B[50] = userPrefp->f302_b; cod_w_list_B[51] = userPrefp->f303_b;
    cod_w_list_B[52] = userPrefp->f310_b; cod_w_list_B[53] = userPrefp->f311_b; cod_w_list_B[54] = userPrefp->f312_b; cod_w_list_B[55] = userPrefp->f313_b;
    cod_w_list_B[56] = userPrefp->f320_b; cod_w_list_B[57] = userPrefp->f321_b; cod_w_list_B[58] = userPrefp->f322_b; cod_w_list_B[59] = userPrefp->f323_b;
    cod_w_list_B[60] = userPrefp->f330_b; cod_w_list_B[61] = userPrefp->f331_b; cod_w_list_B[62] = userPrefp->f332_b; cod_w_list_B[63] = userPrefp->f333_b;

    cod_w_list_C[0] = userPrefp->f000_c; cod_w_list_C[1] = userPrefp->f001_c; cod_w_list_C[2] = userPrefp->f002_c; cod_w_list_C[3] = userPrefp->f003_c;
    cod_w_list_C[4] = userPrefp->f010_c; cod_w_list_C[5] = userPrefp->f011_c; cod_w_list_C[6] = userPrefp->f012_c; cod_w_list_C[7] = userPrefp->f013_c;
    cod_w_list_C[8] = userPrefp->f020_c; cod_w_list_C[9] = userPrefp->f021_c; cod_w_list_C[10] = userPrefp->f022_c; cod_w_list_C[11] = userPrefp->f023_c;
    cod_w_list_C[12] = userPrefp->f030_c; cod_w_list_C[13] = userPrefp->f031_c; cod_w_list_C[14] = userPrefp->f032_c; cod_w_list_C[15] = userPrefp->f033_c;
    cod_w_list_C[16] = userPrefp->f100_c; cod_w_list_C[17] = userPrefp->f101_c; cod_w_list_C[18] = userPrefp->f102_c; cod_w_list_C[19] = userPrefp->f103_c;
    cod_w_list_C[20] = userPrefp->f110_c; cod_w_list_C[21] = userPrefp->f111_c; cod_w_list_C[22] = userPrefp->f112_c; cod_w_list_C[23] = userPrefp->f113_c;
    cod_w_list_C[24] = userPrefp->f120_c; cod_w_list_C[25] = userPrefp->f121_c; cod_w_list_C[26] = userPrefp->f122_c; cod_w_list_C[27] = userPrefp->f123_c;
    cod_w_list_C[28] = userPrefp->f130_c; cod_w_list_C[29] = userPrefp->f131_c; cod_w_list_C[30] = userPrefp->f132_c; cod_w_list_C[31] = userPrefp->f133_c;
    cod_w_list_C[32] = userPrefp->f200_c; cod_w_list_C[33] = userPrefp->f201_c; cod_w_list_C[34] = userPrefp->f202_c; cod_w_list_C[35] = userPrefp->f203_c;
    cod_w_list_C[36] = userPrefp->f210_c; cod_w_list_C[37] = userPrefp->f211_c; cod_w_list_C[38] = userPrefp->f212_c; cod_w_list_C[39] = userPrefp->f213_c;
    cod_w_list_C[40] = userPrefp->f220_c; cod_w_list_C[41] = userPrefp->f221_c; cod_w_list_C[42] = userPrefp->f222_c; cod_w_list_C[43] = userPrefp->f223_c;
    cod_w_list_C[44] = userPrefp->f230_c; cod_w_list_C[45] = userPrefp->f231_c; cod_w_list_C[46] = userPrefp->f232_c; cod_w_list_C[47] = userPrefp->f233_c;
    cod_w_list_C[48] = userPrefp->f300_c; cod_w_list_C[49] = userPrefp->f301_c; cod_w_list_C[50] = userPrefp->f302_c; cod_w_list_C[51] = userPrefp->f303_c;
    cod_w_list_C[52] = userPrefp->f310_c; cod_w_list_C[53] = userPrefp->f311_c; cod_w_list_C[54] = userPrefp->f312_c; cod_w_list_C[55] = userPrefp->f313_c;
    cod_w_list_C[56] = userPrefp->f320_c; cod_w_list_C[57] = userPrefp->f321_c; cod_w_list_C[58] = userPrefp->f322_c; cod_w_list_C[59] = userPrefp->f323_c;
    cod_w_list_C[60] = userPrefp->f330_c; cod_w_list_C[61] = userPrefp->f331_c; cod_w_list_C[62] = userPrefp->f332_c; cod_w_list_C[63] = userPrefp->f333_c;

    
	
	
	return;
}
*/
/*--------------------------------------------------------------------------------*/
/*the main alogorithm involved in the fsim: mutation, recombination, reproduction
 coutfix  are done in that order for the required number of generations
 genNum - the number of generations to be run
 flagReset - flag which says whether to reset the counts for the current call of fsim
 repType - the flag which differentiates between the three types of calls - initrun, prerun and reprun*/
long doFsim(struct userPref_str *userPrefp, long genNum, long flagReset, long repType)
{
    int curSeq, curSite;
    
    getPopScenario(userPrefp, repType, genNum);
    
    for (curGenNum=1; curGenNum <= genNum; curGenNum++)
    {
        printf ("\tGeneration %lu\t\n", curGenNum);
#ifdef DEBUGMODE
        printf("\tGeneration %lu\t\n", curGenNum);
        
        halt;
#endif
        
        
        if (check==0) {
            population_size = userPrefp->prerunSeqNum;
        }
        else if (check==1) {
            if (curGenNum <= userPrefp->uswitch) {
                population_size = userPrefp->initSeqNum;
            }
            else if (curGenNum > userPrefp->uswitch) {
                population_size = userPrefp->SeqNum_b;
            }
            //printf ("Generation = %lu\tpopulation size = %ld\n", curGenNum, population_size);
        }
        /* 150617 matsumoto: record the sequences before mutation. They are used to calculate fitness  */
        for (curSeq=1; curSeq<= population_size; curSeq++) {
            for (curSite=1; curSite <= userPrefp->seqLen; curSite++) {
                prev_seq_dat[curSeq].site[curSite] = seq_dat[curSeq].site[curSite];
            }
        }
        //		curPopSize = popSizeArray[curGenNum];
        doMutations(userPrefp, curGenNum);
        
#ifdef DEBUGMODE
        checkBaseCounts(1);
#endif
        
        doRecombination(userPrefp, curGenNum);
        
#ifdef DEBUGMODE
        checkBaseCounts(2);
#endif
        
        getNextGenSeqs_MCU(userPrefp, curGenNum);			/* HAmod */
        
        updateCounts(userPrefp, curGenNum);
        
#ifdef DEBUGMODE
        checkBaseCounts(3);
#endif
        
        if (repType == INIT_RUN)
        if (curGenNum % 1000 == 0)
        printf("\tGen %ld of %ld done\n", curGenNum, genNum);
    }
    if(flagReset==DO_RESET)
    resetData(userPrefp);
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*Function that fills the popSizeArray for the current cycle. called from
 doFsim at the start of a initrun, prerun or reprun
 */
long getPopScenario(struct userPref_str *userPrefp, long repType, long genNum)
{
    long curPopArrayCell;
    
#ifdef VERBOSE
    printf("\tGetting population scenario\n");
#endif
    switch(repType)
    {
            /*InitRun*/
        case 0:
        {
            for(curPopArrayCell=1;curPopArrayCell<=genNum+1;curPopArrayCell++)
            {
                popSizeArray[curPopArrayCell] = userPrefp->prerunSeqNum;
            }
            break;
        }
            /*PreRun*/
        case 1:
        {
            for(curPopArrayCell=1;curPopArrayCell<=genNum+1;curPopArrayCell++)
            {
                popSizeArray[curPopArrayCell] = userPrefp->initSeqNum;
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
            for(curPopArrayCell=1;curPopArrayCell<=genNum+1;curPopArrayCell++)
            {
                if (curPopArrayCell <= userPrefp->uswitch) {
                    popSizeArray[curPopArrayCell] = userPrefp->initSeqNum;
                }
                else if (curPopArrayCell > userPrefp->uswitch) {
                    popSizeArray[curPopArrayCell] = userPrefp->SeqNum_b;
                }
            }
            break;
        }
    }
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*Function that resets counts and sitedata after every cycle where data is
 not collected. this function is therefore called after initrun and prerun*/
long resetData(struct userPref_str *userPrefp)
{
    long curSite;					/*refers to the current site*/
    
    popFixNum=0;
    popRecNum=0;
    genMultiHitCount=0;
    segMultiHitCount=0;
    for(curSite=1; curSite <= userPrefp->seqLen; curSite++)
    {
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
        siteData[curSite].recNum=0;
    }
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*
 Function that goes through all the sequences finds the expected number of mutations
 and then implements the changes. Uses and updates the polySite array. Resets and uses
 the curSeqHitSite array and curGenHitSite array. curGenHitSite array is reset at the
 start of loop and curSeqHitSite array is reset after it is used. So curSeqHitSite array
 has to be initialized in the init function. Each mutation results in updating the counts
 in the sequence structure .
 
 numbers of segsites and fixations may be a little lower than expected because the number of mutations cannot be higher than the number of sites
 however, this should be encountered so rarely, that it shouldn't have a detectable effect except at very high mutation rates
 
 */
long doMutations(struct userPref_str *userPrefp, long curGenNum)
{
    long curCodPos=0;							/*current codon position, 1, 2=rep, 3=sil*/
    long curSite=0;					/*variable to loop through all flags, also used to implement mutation*/
    long curSeq=0;					/*current sequence being mutated*/
    long curHitSite=0;				/*Index of array storing hitsites*/
    double expHits01=0;             /*expected number of hits for a given codpos*/
    double expHits02=0;
    double expHits03=0;
    double expHits10=0;
    double expHits12=0;
    double expHits13=0;
    double expHits20=0;
    double expHits21=0;
    double expHits23=0;
    double expHits30=0;
    double expHits31=0;
    double expHits32=0;
    
    long numHits01=0;				/*poisson deviate of the expected number of mutations*/
    long numHits02=0;
    long numHits03=0;
    long numHits10=0;
    long numHits12=0;
    long numHits13=0;
    long numHits20=0;
    long numHits21=0;
    long numHits23=0;
    long numHits30=0;
    long numHits31=0;
    long numHits32=0;
    
    long numChanges10=0;				/*effective number of changes from 1 to 0*/
    long numChanges20=0;				/*effective number of changes from 2 to 0*/
    long numChanges30=0;				/*effective number of changes from 3 to 0*/
    long numChanges12=0;				/*effective number of changes from 1 to 2*/
    long numChanges13=0;				/*effective number of changes from 1 to 3*/
    long numChanges23=0;				/*effective number of changes from 2 to 3*/
    long numChanges01=0;				/*effective number of changes from 0 to 1*/
    long numChanges02=0;				/*effective number of changes from 0 to 2*/
    long numChanges03=0;				/*effective number of changes from 0 to 3*/
    long numChanges21=0;				/*effective number of changes from 2 to 1*/
    long numChanges31=0;				/*effective number of changes from 3 to 1*/
    long numChanges32=0;				/*effective number of changes from 3 to 2*/
    
    
    int total_numHits;
    int check2;
    
    //	long i, ct1_pos1, ct1_pos2, ct1_pos3;
    //	double var_d;
    
#ifdef DEBUGMODE
    printf("\tDoing Mutations\n");
    halt;
#endif
    
    /* reset flags for curGenHitSite. curSeqHitsite are reset in getHitSites function - these are malloced in initData() */
    for(curSite=1; curSite <= userPrefp->seqLen; curSite++)
    curGenHitSite[curSite] = 0;
    
    /* Loop over all sequences */
    
    for (curSeq=1; curSeq <= population_size; curSeq++)
    {
        /* Loop over positions 1, 2, 3 of codons */
        for (curCodPos=1; curCodPos<=3; curCodPos++)
        {
            expHits10 = expHits01 = 0.0;
            expHits20 = expHits02 = 0.0;
            expHits30 = expHits03 = 0.0;
            expHits12 = expHits21 = 0.0;
            expHits13 = expHits31 = 0.0;
            expHits23 = expHits32 = 0.0;
            
            /* get expected number of hits for current cod pos */
            /* 150605matsumoto; consider mutations among four nucleotides (total mutation rate * specific mutation rate)*/
            /*if (check==0) {
             check2=0;
             expHits01 = userPrefp->u_a * userPrefp->u01_a * seq_dat[curSeq].baseCount0[curCodPos];
             expHits02 = userPrefp->u_a * userPrefp->u02_a * seq_dat[curSeq].baseCount0[curCodPos];
             expHits03 = userPrefp->u_a * userPrefp->u03_a * seq_dat[curSeq].baseCount0[curCodPos];
             expHits10 = userPrefp->u_a * userPrefp->u10_a * seq_dat[curSeq].baseCount1[curCodPos];
             expHits12 = userPrefp->u_a * userPrefp->u12_a * seq_dat[curSeq].baseCount1[curCodPos];
             expHits13 = userPrefp->u_a * userPrefp->u13_a * seq_dat[curSeq].baseCount1[curCodPos];
             expHits20 = userPrefp->u_a * userPrefp->u20_a * seq_dat[curSeq].baseCount2[curCodPos];
             expHits21 = userPrefp->u_a * userPrefp->u21_a * seq_dat[curSeq].baseCount2[curCodPos];
             expHits23 = userPrefp->u_a * userPrefp->u23_a * seq_dat[curSeq].baseCount2[curCodPos];
             expHits30 = userPrefp->u_a * userPrefp->u30_a * seq_dat[curSeq].baseCount3[curCodPos];
             expHits31 = userPrefp->u_a * userPrefp->u31_a * seq_dat[curSeq].baseCount3[curCodPos];
             expHits32 = userPrefp->u_a * userPrefp->u32_a * seq_dat[curSeq].baseCount3[curCodPos];
             
             numHits10 = numHits01 = 0;
             numHits20 = numHits02 = 0;
             numHits30 = numHits03 = 0;
             numHits12 = numHits21 = 0;
             numHits13 = numHits31 = 0;
             numHits32 = numHits23 = 0;
             }
             else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
             check2=1;
             expHits01 = userPrefp->u_a * userPrefp->u01_a * seq_dat[curSeq].baseCount0[curCodPos];
             expHits02 = userPrefp->u_a * userPrefp->u02_a * seq_dat[curSeq].baseCount0[curCodPos];
             expHits03 = userPrefp->u_a * userPrefp->u03_a * seq_dat[curSeq].baseCount0[curCodPos];
             expHits10 = userPrefp->u_a * userPrefp->u10_a * seq_dat[curSeq].baseCount1[curCodPos];
             expHits12 = userPrefp->u_a * userPrefp->u12_a * seq_dat[curSeq].baseCount1[curCodPos];
             expHits13 = userPrefp->u_a * userPrefp->u13_a * seq_dat[curSeq].baseCount1[curCodPos];
             expHits20 = userPrefp->u_a * userPrefp->u20_a * seq_dat[curSeq].baseCount2[curCodPos];
             expHits21 = userPrefp->u_a * userPrefp->u21_a * seq_dat[curSeq].baseCount2[curCodPos];
             expHits23 = userPrefp->u_a * userPrefp->u23_a * seq_dat[curSeq].baseCount2[curCodPos];
             expHits30 = userPrefp->u_a * userPrefp->u30_a * seq_dat[curSeq].baseCount3[curCodPos];
             expHits31 = userPrefp->u_a * userPrefp->u31_a * seq_dat[curSeq].baseCount3[curCodPos];
             expHits32 = userPrefp->u_a * userPrefp->u32_a * seq_dat[curSeq].baseCount3[curCodPos];
             
             numHits10 = numHits01 = 0;
             numHits20 = numHits02 = 0;
             numHits30 = numHits03 = 0;
             numHits12 = numHits21 = 0;
             numHits13 = numHits31 = 0;
             numHits32 = numHits23 = 0;
             }
             else if (check==1 && (curGenNum > userPrefp->uswitch)) {
             check2=2;
             expHits01 = userPrefp->u_b * userPrefp->u01_b * seq_dat[curSeq].baseCount0[curCodPos];
             expHits02 = userPrefp->u_b * userPrefp->u02_b * seq_dat[curSeq].baseCount0[curCodPos];
             expHits03 = userPrefp->u_b * userPrefp->u03_b * seq_dat[curSeq].baseCount0[curCodPos];
             expHits10 = userPrefp->u_b * userPrefp->u10_b * seq_dat[curSeq].baseCount1[curCodPos];
             expHits12 = userPrefp->u_b * userPrefp->u12_b * seq_dat[curSeq].baseCount1[curCodPos];
             expHits13 = userPrefp->u_b * userPrefp->u13_b * seq_dat[curSeq].baseCount1[curCodPos];
             expHits20 = userPrefp->u_b * userPrefp->u20_b * seq_dat[curSeq].baseCount2[curCodPos];
             expHits21 = userPrefp->u_b * userPrefp->u21_b * seq_dat[curSeq].baseCount2[curCodPos];
             expHits23 = userPrefp->u_b * userPrefp->u23_b * seq_dat[curSeq].baseCount2[curCodPos];
             expHits30 = userPrefp->u_b * userPrefp->u30_b * seq_dat[curSeq].baseCount3[curCodPos];
             expHits31 = userPrefp->u_b * userPrefp->u31_b * seq_dat[curSeq].baseCount3[curCodPos];
             expHits32 = userPrefp->u_b * userPrefp->u32_b * seq_dat[curSeq].baseCount3[curCodPos];
             
             numHits10 = numHits01 = 0;
             numHits20 = numHits02 = 0;
             numHits30 = numHits03 = 0;
             numHits12 = numHits21 = 0;
             numHits13 = numHits31 = 0;
             numHits32 = numHits23 = 0;
             }
             
             
             numHits01 = (long) poidev(expHits01);
             numHits02 = (long) poidev(expHits02);
             numHits03 = (long) poidev(expHits03);
             numHits10 = (long) poidev(expHits10);
             numHits12 = (long) poidev(expHits12);
             numHits13 = (long) poidev(expHits13);
             numHits20 = (long) poidev(expHits20);
             numHits21 = (long) poidev(expHits21);
             numHits23 = (long) poidev(expHits23);
             numHits30 = (long) poidev(expHits30);
             numHits31 = (long) poidev(expHits31);
             numHits32 = (long) poidev(expHits32);
             
             
             #ifdef DEBUGMODE_1
             printf("\n** start loop **\n");
             printf("curSeq: %0ld\tcurCodPos: %0ld\n", curSeq, curCodPos);
             printf("numHits10: %0ld\nnumHits01: %0ld\n\n", numHits10, numHits01);
             printf("numHits20: %0ld\nnumHits02: %0ld\n\n", numHits20, numHits02);
             printf("numHits30: %0ld\nnumHits03: %0ld\n\n", numHits30, numHits03);
             printf("numHits12: %0ld\nnumHits21: %0ld\n\n", numHits12, numHits21);
             printf("numHits13: %0ld\nnumHits31: %0ld\n\n", numHits13, numHits31);
             printf("numHits23: %0ld\nnumHits32: %0ld\n\n", numHits23, numHits32);
             #endif
             
             total_numHits = numHits01 + numHits02 + numHits03 + numHits10 + numHits12 + numHits13 + numHits20 + numHits21 + numHits23 + numHits30 + numHits31 + numHits32;
             
             if (total_numHits > 0)
             {
             if(numHits01 > seq_dat[curSeq].baseCount0[curCodPos]) {numHits01 = seq_dat[curSeq].baseCount0[curCodPos];}
             if(numHits02 > seq_dat[curSeq].baseCount0[curCodPos]) {numHits02 = seq_dat[curSeq].baseCount0[curCodPos];}
             if(numHits03 > seq_dat[curSeq].baseCount0[curCodPos]) {numHits03 = seq_dat[curSeq].baseCount0[curCodPos];}
             if(numHits10 > seq_dat[curSeq].baseCount1[curCodPos]) {numHits10 = seq_dat[curSeq].baseCount1[curCodPos];}
             if(numHits12 > seq_dat[curSeq].baseCount1[curCodPos]) {numHits12 = seq_dat[curSeq].baseCount1[curCodPos];}
             if(numHits13 > seq_dat[curSeq].baseCount1[curCodPos]) {numHits13 = seq_dat[curSeq].baseCount1[curCodPos];}
             if(numHits20 > seq_dat[curSeq].baseCount2[curCodPos]) {numHits20 = seq_dat[curSeq].baseCount2[curCodPos];}
             if(numHits21 > seq_dat[curSeq].baseCount2[curCodPos]) {numHits21 = seq_dat[curSeq].baseCount2[curCodPos];}
             if(numHits23 > seq_dat[curSeq].baseCount2[curCodPos]) {numHits23 = seq_dat[curSeq].baseCount2[curCodPos];}
             if(numHits30 > seq_dat[curSeq].baseCount3[curCodPos]) {numHits30 = seq_dat[curSeq].baseCount3[curCodPos];}
             if(numHits31 > seq_dat[curSeq].baseCount3[curCodPos]) {numHits31 = seq_dat[curSeq].baseCount3[curCodPos];}
             if(numHits32 > seq_dat[curSeq].baseCount3[curCodPos]) {numHits32 = seq_dat[curSeq].baseCount3[curCodPos];}
             */
            
            /* 150612 matsumoto: base count are changed just after the mutation process */
            /* 150613 matsumoto: calculation of numHits should be done just before each mutation (basecounts are changing by each mutation) */
            
            /* mutation from 0 to 1 */
            if (check==0) {
                expHits01 = userPrefp->u_a * userPrefp->u01_a * seq_dat[curSeq].baseCount0[curCodPos];
                numHits01 = 0;
            }
            else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
                expHits01 = userPrefp->u_a * userPrefp->u01_a * seq_dat[curSeq].baseCount0[curCodPos];
                numHits01 = 0;
            }
            else if (check==1 && (curGenNum > userPrefp->uswitch)) {
                expHits01 = userPrefp->u_b * userPrefp->u01_b * seq_dat[curSeq].baseCount0[curCodPos];
                numHits01 = 0;
            }
            /*take poisson deviates of expected numbers */
            numHits01 = (long) poidev(expHits01);
            
            if (numHits01 > 0) {
                if(numHits01 > seq_dat[curSeq].baseCount0[curCodPos]) {
                    numHits01 = seq_dat[curSeq].baseCount0[curCodPos];
                }
                
                getHitSites(userPrefp, '0', seq_dat[curSeq].baseCount0[curCodPos], curSeq, curCodPos, hitSites01, numHits01);
                seq_dat[curSeq].baseCount0[curCodPos] -= numHits01;
                seq_dat[curSeq].baseCount1[curCodPos] += numHits01;
                
                for (curHitSite=1; curHitSite <= numHits01; curHitSite++)
                {
                    
                    curSite = hitSites01[curHitSite];
                    
                    seq_dat[curSeq].site[curSite] = '1';
                    siteData[curSite].mutNum01++;
                    siteData[curSite].freqs0--;
                    siteData[curSite].freqs1++;
                }
            }
            /* mutation from 0 to 2 */
            if (check==0) {
                expHits02 = userPrefp->u_a * userPrefp->u02_a * seq_dat[curSeq].baseCount0[curCodPos];
                numHits02 = 0;
            }
            else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
                expHits02 = userPrefp->u_a * userPrefp->u02_a * seq_dat[curSeq].baseCount0[curCodPos];
                numHits02 = 0;
            }
            else if (check==1 && (curGenNum > userPrefp->uswitch)) {
                expHits02 = userPrefp->u_b * userPrefp->u02_b * seq_dat[curSeq].baseCount0[curCodPos];
                numHits02 = 0;
            }
            /*take poisson deviates of expected numbers */
            numHits02 = (long) poidev(expHits02);
            
            if (numHits02 > 0) {
                if(numHits02 > seq_dat[curSeq].baseCount0[curCodPos]) {
                    numHits02 = seq_dat[curSeq].baseCount0[curCodPos];
                }
                
                
                getHitSites(userPrefp, '0', seq_dat[curSeq].baseCount0[curCodPos], curSeq, curCodPos, hitSites02, numHits02);
                seq_dat[curSeq].baseCount0[curCodPos] -= numHits02;
                seq_dat[curSeq].baseCount2[curCodPos] += numHits02;
                
                for (curHitSite=1; curHitSite <= numHits02; curHitSite++)
                {
                    
                    curSite = hitSites02[curHitSite];
                    
                    seq_dat[curSeq].site[curSite] = '2';
                    siteData[curSite].mutNum02++;
                    siteData[curSite].freqs0--;
                    siteData[curSite].freqs2++;
                }
            }
            /* mutation from 0 to 3 */
            if (check==0) {
                expHits03 = userPrefp->u_a * userPrefp->u03_a * seq_dat[curSeq].baseCount0[curCodPos];
                numHits03 = 0;
            }
            else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
                expHits03 = userPrefp->u_a * userPrefp->u03_a * seq_dat[curSeq].baseCount0[curCodPos];
                numHits03 = 0;
            }
            else if (check==1 && (curGenNum > userPrefp->uswitch)) {
                expHits03 = userPrefp->u_b * userPrefp->u03_b * seq_dat[curSeq].baseCount0[curCodPos];
                numHits03 = 0;
            }
            /*take poisson deviates of expected numbers */
            numHits03 = (long) poidev(expHits03);
            
            if (numHits03 > 0) {
                if(numHits03 > seq_dat[curSeq].baseCount0[curCodPos]) {
                    numHits03 = seq_dat[curSeq].baseCount0[curCodPos];
                }
                
                getHitSites(userPrefp, '0', seq_dat[curSeq].baseCount0[curCodPos], curSeq, curCodPos, hitSites03, numHits03);
                seq_dat[curSeq].baseCount0[curCodPos] -= numHits03;
                seq_dat[curSeq].baseCount3[curCodPos] += numHits03;
                
                for (curHitSite=1; curHitSite <= numHits03; curHitSite++)
                {
                    
                    curSite = hitSites03[curHitSite];
                    
                    seq_dat[curSeq].site[curSite] = '3';
                    siteData[curSite].mutNum03++;
                    siteData[curSite].freqs0--;
                    siteData[curSite].freqs3++;
                }
            }
            /* mutation from 1 to 0 */
            if (check==0) {
                expHits10 = userPrefp->u_a * userPrefp->u10_a * seq_dat[curSeq].baseCount1[curCodPos];
                numHits10 = 0;
            }
            else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
                expHits10 = userPrefp->u_a * userPrefp->u10_a * seq_dat[curSeq].baseCount1[curCodPos];
                numHits10 = 0;
            }
            else if (check==1 && (curGenNum > userPrefp->uswitch)) {
                expHits10 = userPrefp->u_b * userPrefp->u10_b * seq_dat[curSeq].baseCount1[curCodPos];
                numHits10 = 0;
            }
            /*take poisson deviates of expected numbers */
            numHits10 = (long) poidev(expHits10);
            
            if (numHits10 > 0) {
                if(numHits10 > seq_dat[curSeq].baseCount1[curCodPos]) {
                    numHits10 = seq_dat[curSeq].baseCount1[curCodPos];
                }
                
                getHitSites(userPrefp, '1', seq_dat[curSeq].baseCount1[curCodPos], curSeq, curCodPos, hitSites10, numHits10);
                seq_dat[curSeq].baseCount1[curCodPos] -= numHits10;
                seq_dat[curSeq].baseCount0[curCodPos] += numHits10;
                
                for (curHitSite=1; curHitSite <= numHits10; curHitSite++)
                {
                    
                    curSite = hitSites10[curHitSite];
                    
                    seq_dat[curSeq].site[curSite] = '0';
                    siteData[curSite].mutNum10++;
                    siteData[curSite].freqs1--;
                    siteData[curSite].freqs0++;
                }
            }
            /* mutation from 1 to 2 */
            if (check==0) {
                expHits12 = userPrefp->u_a * userPrefp->u12_a * seq_dat[curSeq].baseCount1[curCodPos];
                numHits12 = 0;
            }
            else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
                expHits12 = userPrefp->u_a * userPrefp->u12_a * seq_dat[curSeq].baseCount1[curCodPos];
                numHits12 = 0;
            }
            else if (check==1 && (curGenNum > userPrefp->uswitch)) {
                expHits12 = userPrefp->u_b * userPrefp->u12_b * seq_dat[curSeq].baseCount1[curCodPos];
                numHits12 = 0;
            }
            /*take poisson deviates of expected numbers */
            numHits12 = (long) poidev(expHits12);
            
            if (numHits12 > 0) {
                if(numHits12 > seq_dat[curSeq].baseCount1[curCodPos]) {
                    numHits12 = seq_dat[curSeq].baseCount1[curCodPos];
                }
                
                getHitSites(userPrefp, '1', seq_dat[curSeq].baseCount1[curCodPos], curSeq, curCodPos, hitSites12, numHits12);
                seq_dat[curSeq].baseCount1[curCodPos] -= numHits12;
                seq_dat[curSeq].baseCount2[curCodPos] += numHits12;
                
                for (curHitSite=1; curHitSite <= numHits12; curHitSite++)
                {
                    
                    curSite = hitSites12[curHitSite];
                    
                    seq_dat[curSeq].site[curSite] = '2';
                    siteData[curSite].mutNum12++;
                    siteData[curSite].freqs1--;
                    siteData[curSite].freqs2++;
                }
            }
            /* mutation from 1 to 3 */
            if (check==0) {
                expHits13 = userPrefp->u_a * userPrefp->u13_a * seq_dat[curSeq].baseCount1[curCodPos];
                numHits13 = 0;
            }
            else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
                expHits13 = userPrefp->u_a * userPrefp->u13_a * seq_dat[curSeq].baseCount1[curCodPos];
                numHits13 = 0;
            }
            else if (check==1 && (curGenNum > userPrefp->uswitch)) {
                expHits13 = userPrefp->u_b * userPrefp->u13_b * seq_dat[curSeq].baseCount1[curCodPos];
                numHits13 = 0;
            }
            /*take poisson deviates of expected numbers */
            numHits13 = (long) poidev(expHits13);
            
            if (numHits13 > 0) {
                if(numHits13 > seq_dat[curSeq].baseCount1[curCodPos]) {
                    numHits13 = seq_dat[curSeq].baseCount1[curCodPos];
                }
                
                getHitSites(userPrefp, '1', seq_dat[curSeq].baseCount1[curCodPos], curSeq, curCodPos, hitSites13, numHits13);
                seq_dat[curSeq].baseCount1[curCodPos] -= numHits13;
                seq_dat[curSeq].baseCount3[curCodPos] += numHits13;
                
                for (curHitSite=1; curHitSite <= numHits13; curHitSite++)
                {
                    
                    curSite = hitSites13[curHitSite];
                    
                    seq_dat[curSeq].site[curSite] = '3';
                    siteData[curSite].mutNum13++;
                    siteData[curSite].freqs1--;
                    siteData[curSite].freqs3++;
                }
            }
            /* mutation from 2 to 0 */
            if (check==0) {
                expHits20 = userPrefp->u_a * userPrefp->u20_a * seq_dat[curSeq].baseCount2[curCodPos];
                numHits20 = 0;
            }
            else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
                expHits20 = userPrefp->u_a * userPrefp->u20_a * seq_dat[curSeq].baseCount2[curCodPos];
                numHits20 = 0;
            }
            else if (check==1 && (curGenNum > userPrefp->uswitch)) {
                expHits20 = userPrefp->u_b * userPrefp->u20_b * seq_dat[curSeq].baseCount2[curCodPos];
                numHits20 = 0;
            }
            /*take poisson deviates of expected numbers */
            numHits20 = (long) poidev(expHits20);
            
            if (numHits20 > 0) {
                if(numHits20 > seq_dat[curSeq].baseCount2[curCodPos]) {
                    numHits20 = seq_dat[curSeq].baseCount2[curCodPos];
                }
                
                getHitSites(userPrefp, '2', seq_dat[curSeq].baseCount2[curCodPos], curSeq, curCodPos, hitSites20, numHits20);
                seq_dat[curSeq].baseCount2[curCodPos] -= numHits20;
                seq_dat[curSeq].baseCount0[curCodPos] += numHits20;
                
                for (curHitSite=1; curHitSite <= numHits20; curHitSite++)
                {
                    
                    curSite = hitSites20[curHitSite];
                    
                    seq_dat[curSeq].site[curSite] = '0';
                    siteData[curSite].mutNum20++;
                    siteData[curSite].freqs2--;
                    siteData[curSite].freqs0++;
                }
            }
            /* mutation from 2 to 1 */
            if (check==0) {
                expHits21 = userPrefp->u_a * userPrefp->u21_a * seq_dat[curSeq].baseCount2[curCodPos];
                numHits21 = 0;
            }
            else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
                expHits21 = userPrefp->u_a * userPrefp->u21_a * seq_dat[curSeq].baseCount2[curCodPos];
                numHits21 = 0;
            }
            else if (check==1 && (curGenNum > userPrefp->uswitch)) {
                expHits21 = userPrefp->u_b * userPrefp->u21_b * seq_dat[curSeq].baseCount2[curCodPos];
                numHits21 = 0;
            }
            /*take poisson deviates of expected numbers */
            numHits21 = (long) poidev(expHits21);
            
            if (numHits21 > 0) {
                if(numHits21 > seq_dat[curSeq].baseCount2[curCodPos]) {
                    numHits21 = seq_dat[curSeq].baseCount2[curCodPos];
                }
                
                getHitSites(userPrefp, '2', seq_dat[curSeq].baseCount2[curCodPos], curSeq, curCodPos, hitSites21, numHits21);
                seq_dat[curSeq].baseCount2[curCodPos] -= numHits21;
                seq_dat[curSeq].baseCount1[curCodPos] += numHits21;
                
                for (curHitSite=1; curHitSite <= numHits21; curHitSite++)
                {
                    
                    curSite = hitSites21[curHitSite];
                    
                    seq_dat[curSeq].site[curSite] = '1';
                    siteData[curSite].mutNum21++;
                    siteData[curSite].freqs2--;
                    siteData[curSite].freqs1++;
                }
            }
            /* mutation from 2 to 3 */
            if (check==0) {
                expHits23 = userPrefp->u_a * userPrefp->u23_a * seq_dat[curSeq].baseCount2[curCodPos];
                numHits23 = 0;
            }
            else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
                expHits23 = userPrefp->u_a * userPrefp->u23_a * seq_dat[curSeq].baseCount2[curCodPos];
                numHits23 = 0;
            }
            else if (check==1 && (curGenNum > userPrefp->uswitch)) {
                expHits23 = userPrefp->u_b * userPrefp->u23_b * seq_dat[curSeq].baseCount2[curCodPos];
                numHits23 = 0;
            }
            /*take poisson deviates of expected numbers */
            numHits23 = (long) poidev(expHits23);
            
            if (numHits23 > 0) {
                if(numHits23 > seq_dat[curSeq].baseCount2[curCodPos]) {
                    numHits23 = seq_dat[curSeq].baseCount2[curCodPos];
                }
                
                getHitSites(userPrefp, '2', seq_dat[curSeq].baseCount2[curCodPos], curSeq, curCodPos, hitSites23, numHits23);
                seq_dat[curSeq].baseCount2[curCodPos] -= numHits23;
                seq_dat[curSeq].baseCount3[curCodPos] += numHits23;
                
                for (curHitSite=1; curHitSite <= numHits23; curHitSite++)
                {
                    
                    curSite = hitSites23[curHitSite];
                    
                    seq_dat[curSeq].site[curSite] = '3';
                    siteData[curSite].mutNum23++;
                    siteData[curSite].freqs2--;
                    siteData[curSite].freqs3++;
                }
            }
            /* mutation from 3 to 0 */
            if (check==0) {
                expHits30 = userPrefp->u_a * userPrefp->u30_a * seq_dat[curSeq].baseCount3[curCodPos];
                numHits30 = 0;
            }
            else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
                expHits30 = userPrefp->u_a * userPrefp->u30_a * seq_dat[curSeq].baseCount3[curCodPos];
                numHits30 = 0;
            }
            else if (check==1 && (curGenNum > userPrefp->uswitch)) {
                expHits30 = userPrefp->u_b * userPrefp->u30_b * seq_dat[curSeq].baseCount3[curCodPos];
                numHits30 = 0;
            }
            /*take poisson deviates of expected numbers */
            numHits30 = (long) poidev(expHits30);
            
            if (numHits30 > 0) {
                if(numHits30 > seq_dat[curSeq].baseCount3[curCodPos]) {
                    numHits30 = seq_dat[curSeq].baseCount3[curCodPos];
                }
                
                
                
                getHitSites(userPrefp, '3', seq_dat[curSeq].baseCount3[curCodPos], curSeq, curCodPos, hitSites30, numHits30);
                seq_dat[curSeq].baseCount3[curCodPos] -= numHits30;
                seq_dat[curSeq].baseCount0[curCodPos] += numHits30;
                
                for (curHitSite=1; curHitSite <= numHits30; curHitSite++)
                {
                    
                    curSite = hitSites30[curHitSite];
                    
                    seq_dat[curSeq].site[curSite] = '0';
                    siteData[curSite].mutNum30++;
                    siteData[curSite].freqs3--;
                    siteData[curSite].freqs0++;
                }
            }
            /* mutation from 3 to 1 */
            if (check==0) {
                expHits31 = userPrefp->u_a * userPrefp->u31_a * seq_dat[curSeq].baseCount3[curCodPos];
                numHits31 = 0;
            }
            else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
                expHits31 = userPrefp->u_a * userPrefp->u31_a * seq_dat[curSeq].baseCount3[curCodPos];
                numHits31 = 0;
            }
            else if (check==1 && (curGenNum > userPrefp->uswitch)) {
                expHits31 = userPrefp->u_b * userPrefp->u31_b * seq_dat[curSeq].baseCount3[curCodPos];
                numHits31 = 0;
            }
            /*take poisson deviates of expected numbers */
            numHits31 = (long) poidev(expHits31);
            
            if (numHits31 > 0) {
                if(numHits31 > seq_dat[curSeq].baseCount3[curCodPos]) {
                    numHits31 = seq_dat[curSeq].baseCount3[curCodPos];
                }
                
                
                
                getHitSites(userPrefp, '3', seq_dat[curSeq].baseCount3[curCodPos], curSeq, curCodPos, hitSites31, numHits31);
                seq_dat[curSeq].baseCount3[curCodPos] -= numHits31;
                seq_dat[curSeq].baseCount1[curCodPos] += numHits31;
                
                for (curHitSite=1; curHitSite <= numHits31; curHitSite++)
                {
                    
                    curSite = hitSites31[curHitSite];
                    
                    seq_dat[curSeq].site[curSite] = '1';
                    siteData[curSite].mutNum31++;
                    siteData[curSite].freqs3--;
                    siteData[curSite].freqs1++;
                }
            }
            /* mutation from 3 to 2 */
            if (check==0) {
                expHits32 = userPrefp->u_a * userPrefp->u32_a * seq_dat[curSeq].baseCount3[curCodPos];
                numHits32 = 0;
            }
            else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
                expHits32 = userPrefp->u_a * userPrefp->u32_a * seq_dat[curSeq].baseCount3[curCodPos];
                numHits32 = 0;
            }
            else if (check==1 && (curGenNum > userPrefp->uswitch)) {
                expHits32 = userPrefp->u_b * userPrefp->u32_b * seq_dat[curSeq].baseCount3[curCodPos];
                numHits32 = 0;
            }
            /*take poisson deviates of expected numbers */
            numHits32 = (long) poidev(expHits32);
            
            if (numHits32 > 0) {
                if(numHits32 > seq_dat[curSeq].baseCount3[curCodPos]) {
                    numHits32 = seq_dat[curSeq].baseCount3[curCodPos];
                }
                
                
                
                getHitSites(userPrefp, '3', seq_dat[curSeq].baseCount3[curCodPos], curSeq, curCodPos, hitSites32, numHits32);
                seq_dat[curSeq].baseCount3[curCodPos] -= numHits32;
                seq_dat[curSeq].baseCount2[curCodPos] += numHits32;
                
                for (curHitSite=1; curHitSite <= numHits32; curHitSite++)
                {
                    
                    curSite = hitSites32[curHitSite];
                    
                    seq_dat[curSeq].site[curSite] = '2';
                    siteData[curSite].mutNum32++;
                    siteData[curSite].freqs3--;
                    siteData[curSite].freqs2++;
                }
            }
            
            /*
             if (numHits03 > 0)
             {
             getHitSites(userPrefp, '0', seq_dat[curSeq].baseCount0[curCodPos], curSeq, curCodPos, hitSites03, numHits03);
             seq_dat[curSeq].baseCount0[curCodPos] -= numHits03;
             seq_dat[curSeq].baseCount3[curCodPos] += numHits03;
             
             for (curHitSite=1; curHitSite <= numHits03; curHitSite++)
             {
             
             curSite = hitSites03[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '3';
             siteData[curSite].mutNum03++;
             siteData[curSite].freqs0--;
             siteData[curSite].freqs3++;
             }
             }
             if (numHits12 > 0)
             {
             getHitSites(userPrefp, '1', seq_dat[curSeq].baseCount1[curCodPos], curSeq, curCodPos, hitSites12, numHits12);
             seq_dat[curSeq].baseCount1[curCodPos] -= numHits12;
             seq_dat[curSeq].baseCount2[curCodPos] += numHits12;
             
             for (curHitSite=1; curHitSite <= numHits12; curHitSite++)
             {
             
             curSite = hitSites12[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '2';
             siteData[curSite].mutNum12++;
             siteData[curSite].freqs1--;
             siteData[curSite].freqs2++;
             }
             }
             if (numHits13 > 0)
             {
             getHitSites(userPrefp, '1', seq_dat[curSeq].baseCount1[curCodPos], curSeq, curCodPos, hitSites13, numHits13);
             seq_dat[curSeq].baseCount1[curCodPos] -= numHits13;
             seq_dat[curSeq].baseCount3[curCodPos] += numHits13;
             
             for (curHitSite=1; curHitSite <= numHits13; curHitSite++)
             {
             
             curSite = hitSites13[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '3';
             siteData[curSite].mutNum13++;
             siteData[curSite].freqs1--;
             siteData[curSite].freqs3++;
             }
             }
             if (numHits20 > 0)
             {
             getHitSites(userPrefp, '2', seq_dat[curSeq].baseCount2[curCodPos], curSeq, curCodPos, hitSites20, numHits20);
             seq_dat[curSeq].baseCount2[curCodPos] -= numHits20;
             seq_dat[curSeq].baseCount0[curCodPos] += numHits20;
             
             for (curHitSite=1; curHitSite <= numHits20; curHitSite++)
             {
             
             curSite = hitSites20[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '0';
             siteData[curSite].mutNum20++;
             siteData[curSite].freqs2--;
             siteData[curSite].freqs0++;
             }
             }
             if (numHits21 > 0)
             {
             getHitSites(userPrefp, '2', seq_dat[curSeq].baseCount2[curCodPos], curSeq, curCodPos, hitSites21, numHits21);
             seq_dat[curSeq].baseCount2[curCodPos] -= numHits21;
             seq_dat[curSeq].baseCount1[curCodPos] += numHits21;
             
             for (curHitSite=1; curHitSite <= numHits21; curHitSite++)
             {
             
             curSite = hitSites21[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '1';
             siteData[curSite].mutNum21++;
             siteData[curSite].freqs2--;
             siteData[curSite].freqs1++;
             }
             }
             if (numHits23 > 0)
             {
             getHitSites(userPrefp, '2', seq_dat[curSeq].baseCount2[curCodPos], curSeq, curCodPos, hitSites23, numHits23);
             seq_dat[curSeq].baseCount2[curCodPos] -= numHits23;
             seq_dat[curSeq].baseCount3[curCodPos] += numHits23;
             
             for (curHitSite=1; curHitSite <= numHits23; curHitSite++)
             {
             
             curSite = hitSites23[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '3';
             siteData[curSite].mutNum23++;
             siteData[curSite].freqs2--;
             siteData[curSite].freqs3++;
             }
             }
             if (numHits30 > 0)
             {
             getHitSites(userPrefp, '3', seq_dat[curSeq].baseCount3[curCodPos], curSeq, curCodPos, hitSites30, numHits30);
             seq_dat[curSeq].baseCount3[curCodPos] -= numHits30;
             seq_dat[curSeq].baseCount0[curCodPos] += numHits30;
             
             for (curHitSite=1; curHitSite <= numHits30; curHitSite++)
             {
             
             curSite = hitSites30[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '0';
             siteData[curSite].mutNum30++;
             siteData[curSite].freqs3--;
             siteData[curSite].freqs0++;
             }
             }
             if (numHits31 > 0)
             {
             getHitSites(userPrefp, '3', seq_dat[curSeq].baseCount3[curCodPos], curSeq, curCodPos, hitSites31, numHits31);
             seq_dat[curSeq].baseCount3[curCodPos] -= numHits31;
             seq_dat[curSeq].baseCount1[curCodPos] += numHits31;
             
             for (curHitSite=1; curHitSite <= numHits31; curHitSite++)
             {
             
             curSite = hitSites31[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '1';
             siteData[curSite].mutNum31++;
             siteData[curSite].freqs3--;
             siteData[curSite].freqs1++;
             }
             }
             if (numHits32 > 0)
             {
             getHitSites(userPrefp, '3', seq_dat[curSeq].baseCount3[curCodPos], curSeq, curCodPos, hitSites32, numHits32);
             seq_dat[curSeq].baseCount3[curCodPos] -= numHits32;
             seq_dat[curSeq].baseCount2[curCodPos] += numHits32;
             
             for (curHitSite=1; curHitSite <= numHits32; curHitSite++)
             {
             
             curSite = hitSites32[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '2';
             siteData[curSite].mutNum32++;
             siteData[curSite].freqs3--;
             siteData[curSite].freqs2++;
             }
             }
             
             if (numHits10 > 0)
             {
             #ifdef DEBUGMODE_1
             ct1_pos1 = ct1_pos2 = ct1_pos3 = 0;
             for (i=1; i <= userPrefp->seqLen; i += 3)
             if (seq_dat[curSeq].site[i] == '1')
             ct1_pos1++;
             for (i=2; i <= userPrefp->seqLen; i += 3)
             if (seq_dat[curSeq].site[i] == '1')
             ct1_pos2++;
             for (i=3; i <= userPrefp->seqLen; i += 3)
             if (seq_dat[curSeq].site[i] == '1')
             ct1_pos3++;
             
             printf("numHits10 > 0\n");
             printf("ct1_pos1                 : %0ld\n", ct1_pos1);
             printf("seq_dat[curSeq].baseCount1[1]: %0ld\t%0ld\n", seq_dat[curSeq].baseCount1[1], ct1_pos1 - seq_dat[curSeq].baseCount1[1]);
             printf("ct1_pos2                 : %0ld\n", ct1_pos2);
             printf("seq_dat[curSeq].baseCount1[2]: %0ld\t%0ld\n", seq_dat[curSeq].baseCount1[2], ct1_pos2 - seq_dat[curSeq].baseCount1[2]);
             printf("ct1_pos3                 : %0ld\n", ct1_pos3);
             printf("seq_dat[curSeq].baseCount1[3]: %0ld\t%0ld\n", seq_dat[curSeq].baseCount1[3], ct1_pos3 - seq_dat[curSeq].baseCount1[3]);
             
             printf("10 curSeq=%lu curCodPos=%lu\n", curSeq, curCodPos);
             printf("baseCount1 = %lu baseCount0 = %lu\n", seq_dat[curSeq].baseCount1[curCodPos], seq_dat[curSeq].baseCount0[curCodPos]);
             #endif
             
             getHitSites(userPrefp, '1', seq_dat[curSeq].baseCount1[curCodPos], curSeq, curCodPos, hitSites10, numHits10);
             seq_dat[curSeq].baseCount1[curCodPos] -= numHits10;
             seq_dat[curSeq].baseCount0[curCodPos] += numHits10;
             
             for (curHitSite=1; curHitSite <= numHits10; curHitSite++)
             {
             
             curSite = hitSites10[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '0';
             siteData[curSite].mutNum10++;
             siteData[curSite].freqs1--;
             siteData[curSite].freqs0++;
             }
             
             #ifdef VERBOSE
             printf("10 curSeq=%lu curCodPos=%lu\n", curSeq, curCodPos);
             printf("baseCount1 = %lu baseCount0 = %lu\n", seq_dat[curSeq].baseCount1[curCodPos], seq_dat[curSeq].baseCount0[curCodPos]);
             #endif
             }
             if (numHits01 > 0)
             {
             #ifdef DEBUGMODE_1
             ct1_pos1 = ct1_pos2 = ct1_pos3 = 0;
             for (i=1; i <= userPrefp->seqLen; i += 3)
             if (seq_dat[curSeq].site[i] == '1')
             ct1_pos1++;
             for (i=2; i <= userPrefp->seqLen; i += 3)
             if (seq_dat[curSeq].site[i] == '1')
             ct1_pos2++;
             for (i=3; i <= userPrefp->seqLen; i += 3)
             if (seq_dat[curSeq].site[i] == '1')
             ct1_pos3++;
             
             printf("numHits01 > 0\n");
             printf("ct1_pos1                 : %0ld\n", ct1_pos1);
             printf("seq_dat[curSeq].baseCount1[1]: %0ld\t%0ld\n", seq_dat[curSeq].baseCount1[1], ct1_pos1 - seq_dat[curSeq].baseCount1[1]);
             printf("ct1_pos2                 : %0ld\n", ct1_pos2);
             printf("seq_dat[curSeq].baseCount1[2]: %0ld\t%0ld\n", seq_dat[curSeq].baseCount1[2], ct1_pos2 - seq_dat[curSeq].baseCount1[2]);
             printf("ct1_pos3                 : %0ld\n", ct1_pos3);
             printf("seq_dat[curSeq].baseCount1[3]: %0ld\t%0ld\n", seq_dat[curSeq].baseCount1[3], ct1_pos3 - seq_dat[curSeq].baseCount1[3]);
             
             printf("01 curSeq=%lu curCodPos=%lu\n", curSeq, curCodPos);
             printf("baseCount1 = %lu baseCount0 = %lu\n", seq_dat[curSeq].baseCount1[curCodPos], seq_dat[curSeq].baseCount0[curCodPos]);
             #endif
             
             getHitSites(userPrefp, '0', seq_dat[curSeq].baseCount0[curCodPos], curSeq, curCodPos, hitSites01, numHits01);
             seq_dat[curSeq].baseCount0[curCodPos] -= numHits01;
             seq_dat[curSeq].baseCount1[curCodPos] += numHits01;
             
             for (curHitSite=1; curHitSite <= numHits01; curHitSite++)
             {
             
             curSite = hitSites01[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '1';
             siteData[curSite].mutNum01++;
             siteData[curSite].freqs0--;
             siteData[curSite].freqs1++;
             }
             
             #ifdef VERBOSE
             printf("01 curSeq=%lu curCodPos=%lu\n", curSeq, curCodPos);
             printf("baseCount1 = %lu baseCount0 = %lu\n", seq_dat[curSeq].baseCount1[curCodPos], seq_dat[curSeq].baseCount0[curCodPos]);
             #endif
             }
             
             
             if (numHits02 > 0) {
             for (curHitSite=1; curHitSite <= numHits02; curHitSite++)
             {
             
             curSite = hitSites02[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '2';
             siteData[curSite].mutNum02++;
             siteData[curSite].freqs0--;
             siteData[curSite].freqs2++;
             }
             }
             if (numHits03 > 0) {
             for (curHitSite=1; curHitSite <= numHits03; curHitSite++)
             {
             
             curSite = hitSites03[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '3';
             siteData[curSite].mutNum03++;
             siteData[curSite].freqs0--;
             siteData[curSite].freqs3++;
             }
             }
             if (numHits12 > 0) {
             for (curHitSite=1; curHitSite <= numHits12; curHitSite++)
             {
             
             curSite = hitSites12[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '2';
             siteData[curSite].mutNum12++;
             siteData[curSite].freqs1--;
             siteData[curSite].freqs2++;
             }
             }
             if (numHits13 > 0) {
             for (curHitSite=1; curHitSite <= numHits13; curHitSite++)
             {
             
             curSite = hitSites13[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '3';
             siteData[curSite].mutNum13++;
             siteData[curSite].freqs1--;
             siteData[curSite].freqs3++;
             }
             }
             if (numHits20 > 0) {
             for (curHitSite=1; curHitSite <= numHits20; curHitSite++)
             {
             
             curSite = hitSites20[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '0';
             siteData[curSite].mutNum20++;
             siteData[curSite].freqs2--;
             siteData[curSite].freqs0++;
             }
             }
             if (numHits21 > 0) {
             for (curHitSite=1; curHitSite <= numHits21; curHitSite++)
             {
             
             curSite = hitSites21[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '1';
             siteData[curSite].mutNum21++;
             siteData[curSite].freqs2--;
             siteData[curSite].freqs1++;
             }
             }
             if (numHits23 > 0) {
             for (curHitSite=1; curHitSite <= numHits23; curHitSite++)
             {
             
             curSite = hitSites23[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '3';
             siteData[curSite].mutNum23++;
             siteData[curSite].freqs2--;
             siteData[curSite].freqs3++;
             }
             }
             if (numHits30 > 0) {
             for (curHitSite=1; curHitSite <= numHits30; curHitSite++)
             {
             
             curSite = hitSites30[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '0';
             siteData[curSite].mutNum13++;
             siteData[curSite].freqs3--;
             siteData[curSite].freqs0++;
             }
             }
             if (numHits31 > 0) {
             for (curHitSite=1; curHitSite <= numHits31; curHitSite++)
             {
             
             curSite = hitSites31[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '1';
             siteData[curSite].mutNum31++;
             siteData[curSite].freqs3--;
             siteData[curSite].freqs1++;
             }
             }
             if (numHits32 > 0) {
             for (curHitSite=1; curHitSite <= numHits32; curHitSite++)
             {
             
             curSite = hitSites32[curHitSite];
             
             seq_dat[curSeq].site[curSite] = '2';
             siteData[curSite].mutNum32++;
             siteData[curSite].freqs3--;
             siteData[curSite].freqs2++;
             }
             }
             if (numHits10 > 0) {
             for (curHitSite=1; curHitSite <= numHits10; curHitSite++)
             {
             #ifdef DEBUGMODE_1
             printf("$$ numHits10 loop\n");
             #endif
             
             curSite = hitSites10[curHitSite];
             
             #ifdef VERBOSE
             printf("%lu\n", curSite);
             #endif
             #ifdef DEBUGMODE
             if (seq_dat[curSeq].site[curSite] != '1')
             printf("10 not 1\n");
             #endif
             
             #ifdef DEBUGMODE_1
             if (seq_dat[curSeq].site[curSite] != '1')
             {
             printf("mutation error 10 hit at -%c- site\n", seq_dat[curSeq].site[curSite] );
             halt;
             }
             #endif
             
             seq_dat[curSeq].site[curSite] = '0';
             siteData[curSite].mutNum10++;
             siteData[curSite].freqs1--;
             siteData[curSite].freqs0++;
             
             
             //	sArray[curSite]=sampleSval(curSite, 1);
             
             }
             }
             
             if (numHits01 > 0) {
             for (curHitSite=1; curHitSite <= numHits01; curHitSite++)
             {
             //printf("$$ numHits01 loop\n");
             
             curSite = hitSites01[curHitSite];
             
             #ifdef VERBOSE
             printf("%lu\n", curSite);
             #endif
             #ifdef DEBUGMODE
             if (seq_dat[curSeq].site[curSite] !='0')
             printf("01 not 0\n");
             #endif
             
             #ifdef DEBUGMODE_1
             if (seq_dat[curSeq].site[curSite] != '0')
             {
             printf("mutation error 01 hit at -%c- site\n", seq_dat[curSeq].site[curSite] );
             halt;
             }
             #endif
             
             seq_dat[curSeq].site[curSite] = '1';
             siteData[curSite].mutNum01++;
             siteData[curSite].freqs0--;
             siteData[curSite].freqs1++;
             
             //  sArray[curSite]=sampleSval(curSite, 0);
             
             }
             }*/
            /* 150612 matsumoto: for more tha two alleles, base count change should be done just after the getHitSite process */
            
            //	printf("!!!  processing basect changes !!!\n");
            /*numChanges10 = 0;
             numChanges10 = numHits10 - numHits01;
             seq_dat[curSeq].baseCount1[curCodPos] -= numChanges10;
             seq_dat[curSeq].baseCount0[curCodPos] += numChanges10;
             
             numChanges20 = 0;
             numChanges20 = numHits20 - numHits02;
             seq_dat[curSeq].baseCount2[curCodPos] -= numChanges20;
             seq_dat[curSeq].baseCount0[curCodPos] += numChanges20;
             
             numChanges30 = 0;
             numChanges30 = numHits30 - numHits03;
             seq_dat[curSeq].baseCount3[curCodPos] -= numChanges30;
             seq_dat[curSeq].baseCount0[curCodPos] += numChanges30;
             
             numChanges12 = 0;
             numChanges12 = numHits12 - numHits21;
             seq_dat[curSeq].baseCount1[curCodPos] -= numChanges12;
             seq_dat[curSeq].baseCount2[curCodPos] += numChanges21;
             
             numChanges13 = 0;
             numChanges13 = numHits13 - numHits31;
             seq_dat[curSeq].baseCount1[curCodPos] -= numChanges13;
             seq_dat[curSeq].baseCount3[curCodPos] += numChanges31;
             
             numChanges23 = 0;
             numChanges23 = numHits23 - numHits32;
             seq_dat[curSeq].baseCount2[curCodPos] -= numChanges23;
             seq_dat[curSeq].baseCount3[curCodPos] += numChanges32;*/
            
            
            
            /*in DEBUGMODE, output the information of only 1->0 and 0->1 changes*/
#ifdef VERBOSE
            if(numChanges10)
            printf("baseCount1 = %lu baseCount0 = %lu\n", seq_dat[curSeq].baseCount1[curCodPos], seq_dat[curSeq].baseCount0[curCodPos]);
#endif
            
#ifdef DEBUGMODE_1
            ct1_pos1 = ct1_pos2 = ct1_pos3 = 0;
            for (i=1; i <= userPrefp->seqLen; i += 3)
            if (seq_dat[curSeq].site[i] == '1')
            ct1_pos1++;
            for (i=2; i <= userPrefp->seqLen; i += 3)
            if (seq_dat[curSeq].site[i] == '1')
            ct1_pos2++;
            for (i=3; i <= userPrefp->seqLen; i += 3)
            if (seq_dat[curSeq].site[i] == '1')
            ct1_pos3++;
            
            printf("post doMutations update:\n");
            printf("curSeq: %0ld\tcurCodPos: %0ld\n", curSeq, curCodPos);
            printf("numHits10: %0ld\nnumHits01: %0ld\n", numHits10, numHits01);
            if ((numHits10 > 0) || (numHits01 > 0))
            printf("numChanges10: %0ld\n", numChanges10);
            
            printf("ct1_pos1                 : %0ld\n", ct1_pos1);
            printf("seq_dat[curSeq].baseCount1[1]: %0ld\t%0ld\t", seq_dat[curSeq].baseCount1[1], ct1_pos1 - seq_dat[curSeq].baseCount1[1]);
            if (ct1_pos1 != seq_dat[curSeq].baseCount1[1])
            printf("******\n");
            else
            printf("\n");
            printf("ct1_pos2                 : %0ld\n", ct1_pos2);
            printf("seq_dat[curSeq].baseCount1[2]: %0ld\t%0ld\t", seq_dat[curSeq].baseCount1[2], ct1_pos2 - seq_dat[curSeq].baseCount1[2]);
            if (ct1_pos2 != seq_dat[curSeq].baseCount1[2])
            printf("******\n");
            else
            printf("\n");
            printf("ct1_pos3                 : %0ld\n", ct1_pos3);
            printf("seq_dat[curSeq].baseCount1[3]: %0ld\t%0ld\t", seq_dat[curSeq].baseCount1[3], ct1_pos3 - seq_dat[curSeq].baseCount1[3]);
            if (ct1_pos3 != seq_dat[curSeq].baseCount1[3])
            printf("******\n");
            else
            printf("\n");
#endif
        }
    }
    
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*
 Searches for the sites that could be mutated. This is where the
 multihit restrictions are implemented. If nomultihit flag is on, then no multiple hits
 within a generation and at segregating sites are permitted. Regardless of nomultihit,
 the same site in a given gene cannot be hit twice in a given generation.
 The function also checks every 'limit' passes through the loop if any more mutations are possible.
 This is to prevent cases where the flags prevent all mutations to occur
 
 Parameters passed
 origSite 	- The kind of base that gets mutated (1 or 0)
 codPos		- The codon position 1, 2 repsites 3 silsite
 baseCount	- Total number of such bases at all codPos sites
 curSeq		- The current sequence where the search is made
 hitSites	- The array which is filled with positions where mutations should occur
 numHits	 	- The number of hits that should occur
 
 HA - important that when multiple hits are encountered, the code will search for a new site for the mutation
 numbers of segsites and fixations should not be lowered due to how multiple hits are processed
 
 no multiple mutations at the same site in the same chr in the same generation (independent of userPrefp->noMultiHit)
 userPrefp->noMultiHit == 1
	no multiple mutations at segregating sites
 */
long getHitSites(struct userPref_str *userPrefp, char origSite, long baseCount, long curSeq, long codPos, long* hitSites, long numHits)
{
    long limit;					/*The maximum number of iterations through the loop
                                 before checking if there are actually any sites left
                                 that can be mutated. If not then error*/
    long numRandSearch;			/*the number of searches made for the current function call*/
    long curNumHits;				/*number of hits located at any time*/
    long hitCod;					/*a mutation is possible at the hitCod codon with
                                     origSite at codPos*/
    long curSite;					/*The actual base position for the hit*/
    long numUnflagged, numFlagged;				/*The number of unflagged */
    
    /*limit is the number of random searches made before a check is made on
     the sites to see if the remaining mutations are posssible*/
    limit = numHits * MULTIHITFACTOR;
    numRandSearch=0;
    curNumHits=0;
    
    /*repeat until all hitsites are found*/
    while (curNumHits < numHits)
    {
        /*keep track of the number of random searches*/
        numRandSearch++;
        
        /*if the number of searches is greater than or equal to a set limit*/
        if (numRandSearch >= limit)
        {
            /*count the number of unflagged sites*/
            numUnflagged = numFlagged = 0;
            if (userPrefp->noMultiHit == 1)
            {
                /* if noMultihit is on then check all sites if the polysites flag
                 is set. this is because curGenHitSites are also flagged as polySite */
                for (curSite = codPos; curSite <= userPrefp->seqLen; curSite+=3)
                if (seq_dat[curSeq].site[curSite] == origSite)
                {
                    if (polySite[curSite]==0)
                    numUnflagged++;
                    else
                    numFlagged++;
                }
            }
            else
            {
                /* if noMultiHit is off then check all sites if curSeqHitSite flag is set */
                for (curSite=codPos; curSite <= userPrefp->seqLen; curSite+=3)
                if (seq_dat[curSeq].site[curSite] == origSite)
                {
                    if (curSeqHitSite[curSite] == 0)
                    numUnflagged++;
                    else
                    numFlagged++;
                }
            }
            if (numUnflagged < (numHits - curNumHits))
            {
                for (curSite = codPos; curSite <= userPrefp->seqLen; curSite+=3) {
                    printf ("%c", seq_dat[curSeq].site[curSite]);
                }
                printf ("\n");
                printf("Very High mutation rate / site\n");
                printf("origSite    : %c\n", origSite);
                printf("codPos      : %0ld\n", codPos);
                printf("basecount0  : %0ld\n", seq_dat[curSeq].baseCount0[codPos]);
                //printf("tot_1_ct_silent    : %0ld\n", seq_dat[curSeq].tot_1_ct_silent);
                //printf("tot_1_ct_replacement    : %0ld\n", seq_dat[curSeq].tot_1_ct_replacement);
                printf("numFlagged  : %0ld\n", numFlagged);
                printf("numUnflagged: %0ld\n", numUnflagged);
                printf("numHits     : %0ld\n", numHits);
                printf("curNumHits  : %0ld\n", curNumHits);
                printf("please reduce mutation rate or remove nomultihit flag\n");
                exit(1);
            }
            else
            {
                /* never repeat this process again for the current sequence set limit to maximum value */
                limit = ULONG_MAX;
            }
        }
        /* find the site where a mutation can occur */
        hitCod = getRandLong_(1, baseCount);
        curSite = getPosNthChar_v2(userPrefp, seq_dat[curSeq].site, curSeq, codPos, hitCod, origSite);
        /*
         seq_dat[curSeq].site
         &(seq_dat[curSeq].site)
         &(seq_dat[curSeq].site[0])
         */
        
#ifdef DEBUGMODE
        if(seq_dat[curSeq].site[curSite] != origSite)
        printf("getPosNthChar returned site which is not correct");
#endif
        
        /*if multiple hits are not permitted then check multiple hit rules and if violated go back and get another site*/
        if (userPrefp->noMultiHit == 1)
        {
            /*The order in which these checks are made is significant because
             curGenHitSites are always polySites but the reverse is not true*/
            if (curGenHitSite[curSite] == 1)
            {
                genMultiHitCount++;
                continue;
            }
            else if (polySite[curSite] == 1)
            {
                segMultiHitCount++;
                continue;
            }
        }
        /* No matter the settings for multiple hit preferences, multiple hits in the same sequence at same site are not allowed */
        if (curSeqHitSite[curSite] == 0)
        {
            /*update counts for the while loop*/
            curNumHits++;
            
            /*set flags for polySite, curSeqHit and curGenHit*/
            polySite[curSite] = 1;
            curSeqHitSite[curSite] = 1;
            curGenHitSite[curSite] = 1;
            
            /*store the curSite in the hitSites array*/
            hitSites[curNumHits] = curSite;
        }
    }
    
    /*resets the curSeqHitSite array. Since it is reset after it is used first
     the array has to be initialized with all zeros in the initdata function*/
    for (curSite=1; curSite <= numHits; curSite++)
    {
        curSeqHitSite[hitSites[curSite]] = 0;
    }
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*
 searches for the nth 1 or 0 for a given codon position. for example if
 codpos=2, origsite=1 and hitCod=20, then the function will return the site position of
 the 20th codon with a 1 at its position 2
 site			- the character array storing the given sequence
 codPos		- the position in the codon ie 1, 2 or 3
 origSite		- the function will search for codons with origSite at codPos
 hitCod		- the codon which has been hit
 
 HA origSite was treated as int/long although it's passed as char - changed to char comparisons/outputs
 */
long getPosNthChar(struct userPref_str *userPrefp, char *site,  long codPos, long hitCod, char origSite)
{
    long curCod;					/*refers to current codon*/
    long curSite, found;					/*refers to current site*/
    
    
    curCod = 0;
    curSite = codPos;
    /*loop until curCod is the hitCod codon with origSite at the codPos
     The check for curSite is redundant because we already know that there
     are enough codons to satisfy the condition. So this should just serve
     as an error check. The second condition may be removed */
    found = 0;
    while ((curSite <= userPrefp->seqLen) && (found == 0))
    {
        if (site[curSite] == origSite)
        {
            curCod++;
            if (curCod == hitCod)
            found = 1;
        }
        /*increment current site by 3 to look at codPos in the next codon*/
        if (found == 0)
        curSite += 3;
    }
    
    
    /*The error check. This may be removed */			/* HA_speed_2 */
#ifdef DEBUGMODE_1
    if (found == 0)
    {
        printf("getPosNthChar\n");
        printf("\tsite all:\n");
        for (curSite=1; curSite <= userPrefp->seqLen; curSite++)
        printf("%c", site[curSite]);
        printf("-end-\n");
        
        printf("\tsite pos 1: \n");
        for (curSite=1; curSite <= userPrefp->seqLen; curSite += 3)
        printf("%c\t%0ld\n", site[curSite], curSite);
        printf("-end-\n");
        
        printf("\tsite pos 2: \n");
        for (curSite=2; curSite <= userPrefp->seqLen; curSite += 3)
        printf("%c\t%0ld\n", site[curSite], curSite);
        printf("-end-\n");
        
        printf("\tsite pos 3: \n");
        for (curSite=3; curSite <= userPrefp->seqLen; curSite += 3)
        printf("%c\t%0ld\n", site[curSite], curSite);
        printf("-end-\n");
        
        printf("\tcodPos  : %0ld\n", codPos);
        printf("\thitCod  : %0ld\n", hitCod);
        printf("\torigSite: %c\n", origSite);
        
        printf("looking for %ldth %c\n", hitCod, origSite);
        curCod=0;
        for (curSite=codPos; curSite <= userPrefp->seqLen; curSite += 3)
        {
            if (site[curSite] == origSite)
            curCod++;
            printf("%c\t%ld\t%ld\n", site[curSite], curSite, curCod);
        }
        printf("\nfound a total of %lu %c's in the sequence\n", curCod, origSite);
        errorOut(("Something wrong with baseCounts. cursite>seqlen. Check Code\n"));
    }
#endif
    
    //	printf("\tcurSite: %ld\n", curSite);
    return curSite;
}
/*--------------------------------------------------------------------------------*/
long getPosNthChar_v2(struct userPref_str *userPrefp, char *site, long curSeq, long codPos, long hitCod, char origSite)
{
    long curCod;					/*refers to current codon*/
    long curSite;					/*refers to current site*/
    
    curCod = 0;
    curSite = codPos;
    /*loop until curCod is the hitCod codon with origSite at the codPos
     The check for curSite is redundant because we already know that there
     are enough codons to satisfy the condition. So this should just serve
     as an error check. The second condition may be removed */
    while (curSite <= userPrefp->seqLen)
    {
        //printf("%ld\n", codPos);
        //printf("%c,", seq_dat[curSeq].site[curSite]);
        if (site[curSite] == origSite)
        {
            //printf ("%ld, %ld\n", curCod, curSite);
            curCod++;
            if (curCod == hitCod)
            return curSite;
        }
        /*increment current site by 3 to look at codPos in the next codon*/
        curSite += 3;
    }
    printf("\ngetPosNthChar_v2 error\n");
    printf("\tfound a total of %lu %c's in the sequence: looking for %0ldth occurrence %ld\n", curCod, origSite, hitCod, codPos);
    errorOut(("\tSomething wrong with baseCounts. cursite>seqlen. Check Code\n"));
    
    return curSite;
}
/*--------------------------------------------------------------------------------*/
/*Function that goes through all the sequences, randomly pairing them and performs
 recombination between the paired sequences. Values are changed only when there is
 actual difference between recombining sites.
 the following are malloced in initData()
 long seqIndex[MAXSEQNUM+1];
 long crossOverSite[MAXSEQLEN+1];
 */
long doRecombination(struct userPref_str *userPrefp, long curGenNum)
{
    double expCrossOvers;					/*The expected number of crossovers*/
    long maxCrossOvers;			/*The maximum number of crossovers = seqLen-1*/
    long numCrossOvers;			/*The number of crossovers for the current pair*/
    long curSeq;					/*Index of the first of the recombining seq in the seqIndex array*/
    long curSeqPlusOne;			/*Index of the second recombining seq in the seqIndex array*/
    long curRecSeq1;				/*Sequence number of the current recombining sequence*/
    long curRecSeq2;				/*Sequence number of the current recombining sequence*/
    long nextSeq;					/*Variable used to locate the second parent randombly*/
    long tempL;					/*temporary long variable*/
    long curCrossSiteNum;			/*The number of cross over sites identified so far*/
    long curCrossSite;				/*currently identified site for crossover*/
    long curSeg;					/*The current segment which is being exchanged*/
    long segStart;					/*Site at which the current segment starts*/
    long segEnd;					/*Site before which the current segment ends*/
    long curSite;					/*The current site which is being exchanged*/
    long curCodPos;							/*codPos of the current site*/
    //	long i, j, baseCount0, baseCount1;
    
#ifdef DEBUGMODE
    printf("\tDoing Recombination\n");
    halt;
#endif
    
    
    /*expected number of crossovers is the same for all pairs of sequences */
    if (check==0) {
        population_size = userPrefp->prerunSeqNum;
        expCrossOvers = userPrefp->c * userPrefp->seqLen;
    }
    else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
        population_size = userPrefp->initSeqNum;
        expCrossOvers = userPrefp->c * userPrefp->seqLen;
        //printf ("recombination = %f\n", expCrossOvers);
    }
    else if (check==1 && (curGenNum > userPrefp->uswitch)) {
        population_size = userPrefp->SeqNum_b;
        expCrossOvers = userPrefp->c_b * userPrefp->seqLen;
        //printf ("recombination = %f\n", expCrossOvers);
    }
    
    
    /*The maximum number of crossovers is one less than the sequence length
     so if there are only 2 bases in the sequence the maximum number of
     crossovers =1*/
    maxCrossOvers = userPrefp->seqLen - 1;
    
    /*seqIndex array is used to pair up sequences. It is initialiezed with
     numbers 1 to userPrefp->initSeqNum*/
    
    for(curSeq=1; curSeq<= population_size; curSeq++)
    seqIndex[curSeq] = curSeq;
    
    /*The indices of the pairs of sequences that are going to be recombined are
     stored at positions 1 and 2 in the seqIndex array*/
    curSeq = 1;
    curSeqPlusOne = 2;
    
    /*repeat loop until the second of the pair is the last cell in the array*/
    while (curSeqPlusOne <= population_size)
    {
        /*find an index to be swapped with the index at curSeqPlusOne position in the
         seqIndex array. Randomly pick a cell from all the remaining indices*/
        nextSeq = getRandLong_(curSeqPlusOne, population_size);
        /*swap the index from the newly located site with that at curSeqPlusOne*/
        tempL = seqIndex[curSeqPlusOne];
        seqIndex[curSeqPlusOne] = seqIndex[nextSeq];
        seqIndex[nextSeq] = tempL;
        /*The indices of the random pair of sequences are at curSeq and curSeqPlusOne*/
        curRecSeq1 = seqIndex[curSeq];
        curRecSeq2 = seqIndex[curSeqPlusOne];
#ifdef VERBOSE
        printf("%lu-recombing %lu and %lu\n", curSeq, curRecSeq1, curRecSeq2);
        halt;
#endif
        /*update curSeq and curSeqPlusOne to point to the next pair of cells in the seqIndex array*/
        curSeq += 2;
        curSeqPlusOne += 2;
        /*The number of crossovers is the poisson deviate limited by the maxCrossOvers*/
        numCrossOvers = poidev(expCrossOvers);
        if (numCrossOvers < 0)
        numCrossOvers = 0;
        if(numCrossOvers > maxCrossOvers)
        numCrossOvers = maxCrossOvers;
        popRecNum += numCrossOvers;
        /*If there is going to be any crossOvers then process*/
        if (numCrossOvers != 0)
        {
            /*crossOverSite corresponds to the site that marks the beginning of
             each segement. so the first segment always starts at 1*/
            crossOverSite[1]=1;
            /*generate numCrossOvers sites using crossHitSite array of flags to
             identify multiple hits. store the sites in crossOverSite array*/
            curCrossSiteNum=1;
            while (curCrossSiteNum<=numCrossOvers)
            {
                curCrossSite=getRandLong_(2, userPrefp->seqLen);
                if (crossHitSite[curCrossSite] == 0)
                {
                    /*curCrossSiteNum is incremented first before storing the value because
                     even when there is only one crossover we take it as 1+1+1 to consider
                     the beginning and end of the sequence as two crossover sites. So we will
                     have crooOverSite=[1, crossSite, seqLen].numCrossOvers will be incremented
                     by 2 after flags are reset
                     */
                    curCrossSiteNum++;
                    crossOverSite[curCrossSiteNum]=curCrossSite;
                    crossHitSite[curCrossSite] = 1;
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
            for(curSeg=1; curSeg<=userPrefp->seqLen; curSeg++)
            {
                if(crossHitSite[curSeg]==1)
                {
                    printf("%lu ", curSeg);
                    getchar();getchar();
                }
            }
            /*numCrossOvers is incremented by 2 to account for the beginning and end*/
            numCrossOvers = numCrossOvers + 1;
            crossOverSite[numCrossOvers] = userPrefp->seqLen + 1;
            /*sort the crossOverSite array. first and last positions are already sorted
             so we are only sorting the others. so sort only if more than 1 site other than
             the beginning and end ie only if numCrossOvers>3*/
            if(numCrossOvers>3)
            qsort(&crossOverSite[2], numCrossOvers-1, sizeof(long), compareUL);
#ifdef VERBOSE
            printf("Found all crossover sites\n");
            halt;
            for(curSeg=1; curSeg<=numCrossOvers; curSeg++)
            printf("%lu\t", crossOverSite[curSeg]);
            ent;
#endif
            /*Go through alternate segments starting from segment 1. Loop until
             curSeg>numCrossOvers-1*/
            curSeg=1;
            while (curSeg < numCrossOvers)
            {
                /*A segment starts from one crossOverSite to next crossOverSite-1
                 this is why the <segEng check in for loop*/
                segStart = crossOverSite[curSeg];
                segEnd = crossOverSite[curSeg+1];
#ifdef VERBOSE
                printf("Exchanging segment %lu-%lu\n", segStart, segEnd-1);
#endif
                for(curSite = segStart; curSite < segEnd; curSite++)
                {
                    /*Even for segments that get exchanged, swap sites only if
                     they are different*/
                    if (polySite[curSite])
                    {
                        if(seq_dat[curRecSeq1].site[curSite] == '0')
                        {
                            if (seq_dat[curRecSeq2].site[curSite] == '1')
                            {
                                /*site in first seq is 0 and in second is 1
                                 so swap and update counts*/
                                curCodPos=sitePosToCodPos(curSite);
                                seq_dat[curRecSeq1].site[curSite] ='1';
                                seq_dat[curRecSeq1].baseCount0[curCodPos]--;
                                seq_dat[curRecSeq1].baseCount1[curCodPos]++;
                                seq_dat[curRecSeq2].site[curSite] = '0';
                                seq_dat[curRecSeq2].baseCount0[curCodPos]++;
                                seq_dat[curRecSeq2].baseCount1[curCodPos]--;
                            }
                            else if (seq_dat[curRecSeq2].site[curSite] == '2')
                            {
                                /*site in first seq is 0 and in second is 2
                                 so swap and update counts*/
                                curCodPos=sitePosToCodPos(curSite);
                                seq_dat[curRecSeq1].site[curSite] ='2';
                                seq_dat[curRecSeq1].baseCount0[curCodPos]--;
                                seq_dat[curRecSeq1].baseCount2[curCodPos]++;
                                seq_dat[curRecSeq2].site[curSite] = '0';
                                seq_dat[curRecSeq2].baseCount0[curCodPos]++;
                                seq_dat[curRecSeq2].baseCount2[curCodPos]--;
                            }
                            else if (seq_dat[curRecSeq2].site[curSite] == '3')
                            {
                                /*site in first seq is 0 and in second is 3
                                 so swap and update counts*/
                                curCodPos=sitePosToCodPos(curSite);
                                seq_dat[curRecSeq1].site[curSite] = '3';
                                seq_dat[curRecSeq1].baseCount0[curCodPos]--;
                                seq_dat[curRecSeq1].baseCount3[curCodPos]++;
                                seq_dat[curRecSeq2].site[curSite] = '0';
                                seq_dat[curRecSeq2].baseCount0[curCodPos]++;
                                seq_dat[curRecSeq2].baseCount3[curCodPos]--;
                            }
                        }
                        
                        else if(seq_dat[curRecSeq1].site[curSite] == '1')
                        {
                            if (seq_dat[curRecSeq2].site[curSite] == '0')
                            {
                                /*site in first seq is 1 and in second is 0
                                 so swap and update counts*/
                                curCodPos=sitePosToCodPos(curSite);
                                seq_dat[curRecSeq1].site[curSite] ='0';
                                seq_dat[curRecSeq1].baseCount1[curCodPos]--;
                                seq_dat[curRecSeq1].baseCount0[curCodPos]++;
                                seq_dat[curRecSeq2].site[curSite] = '1';
                                seq_dat[curRecSeq2].baseCount1[curCodPos]++;
                                seq_dat[curRecSeq2].baseCount0[curCodPos]--;
                            }
                            else if (seq_dat[curRecSeq2].site[curSite] == '2')
                            {
                                /*site in first seq is 1 and in second is 2
                                 so swap and update counts*/
                                curCodPos=sitePosToCodPos(curSite);
                                seq_dat[curRecSeq1].site[curSite] ='2';
                                seq_dat[curRecSeq1].baseCount1[curCodPos]--;
                                seq_dat[curRecSeq1].baseCount2[curCodPos]++;
                                seq_dat[curRecSeq2].site[curSite] = '1';
                                seq_dat[curRecSeq2].baseCount1[curCodPos]++;
                                seq_dat[curRecSeq2].baseCount2[curCodPos]--;
                            }
                            else if (seq_dat[curRecSeq2].site[curSite] == '3')
                            {
                                /*site in first seq is 1 and in second is 3
                                 so swap and update counts*/
                                curCodPos=sitePosToCodPos(curSite);
                                seq_dat[curRecSeq1].site[curSite] ='3';
                                seq_dat[curRecSeq1].baseCount1[curCodPos]--;
                                seq_dat[curRecSeq1].baseCount3[curCodPos]++;
                                seq_dat[curRecSeq2].site[curSite] = '1';
                                seq_dat[curRecSeq2].baseCount1[curCodPos]++;
                                seq_dat[curRecSeq2].baseCount3[curCodPos]--;
                            }
                        }
                        
                        else if(seq_dat[curRecSeq1].site[curSite] == '2')
                        {
                            if (seq_dat[curRecSeq2].site[curSite] == '0')
                            {
                                /*site in first seq is 2 and in second is 0
                                 so swap and update counts*/
                                curCodPos=sitePosToCodPos(curSite);
                                seq_dat[curRecSeq1].site[curSite] ='0';
                                seq_dat[curRecSeq1].baseCount2[curCodPos]--;
                                seq_dat[curRecSeq1].baseCount0[curCodPos]++;
                                seq_dat[curRecSeq2].site[curSite] = '2';
                                seq_dat[curRecSeq2].baseCount2[curCodPos]++;
                                seq_dat[curRecSeq2].baseCount0[curCodPos]--;
                            }
                            else if (seq_dat[curRecSeq2].site[curSite] == '1')
                            {
                                /*site in first seq is 2 and in second is 1
                                 so swap and update counts*/
                                curCodPos=sitePosToCodPos(curSite);
                                seq_dat[curRecSeq1].site[curSite] ='1';
                                seq_dat[curRecSeq1].baseCount2[curCodPos]--;
                                seq_dat[curRecSeq1].baseCount1[curCodPos]++;
                                seq_dat[curRecSeq2].site[curSite] = '2';
                                seq_dat[curRecSeq2].baseCount2[curCodPos]++;
                                seq_dat[curRecSeq2].baseCount1[curCodPos]--;
                            }
                            else if (seq_dat[curRecSeq2].site[curSite] == '3')
                            {
                                /*site in first seq is 2 and in second is 3
                                 so swap and update counts*/
                                curCodPos=sitePosToCodPos(curSite);
                                seq_dat[curRecSeq1].site[curSite] ='3';
                                seq_dat[curRecSeq1].baseCount2[curCodPos]--;
                                seq_dat[curRecSeq1].baseCount3[curCodPos]++;
                                seq_dat[curRecSeq2].site[curSite] = '2';
                                seq_dat[curRecSeq2].baseCount2[curCodPos]++;
                                seq_dat[curRecSeq2].baseCount3[curCodPos]--;
                            }
                        }
                        
                        else if(seq_dat[curRecSeq1].site[curSite] == '3')
                        {
                            if (seq_dat[curRecSeq2].site[curSite] == '0')
                            {
                                /*site in first seq is 3 and in second is 0
                                 so swap and update counts*/
                                curCodPos=sitePosToCodPos(curSite);
                                seq_dat[curRecSeq1].site[curSite] ='0';
                                seq_dat[curRecSeq1].baseCount3[curCodPos]--;
                                seq_dat[curRecSeq1].baseCount0[curCodPos]++;
                                seq_dat[curRecSeq2].site[curSite] = '3';
                                seq_dat[curRecSeq2].baseCount3[curCodPos]++;
                                seq_dat[curRecSeq2].baseCount0[curCodPos]--;
                            }
                            else if (seq_dat[curRecSeq2].site[curSite] == '1')
                            {
                                /*site in first seq is 3 and in second is 1
                                 so swap and update counts*/
                                curCodPos=sitePosToCodPos(curSite);
                                seq_dat[curRecSeq1].site[curSite] ='1';
                                seq_dat[curRecSeq1].baseCount3[curCodPos]--;
                                seq_dat[curRecSeq1].baseCount1[curCodPos]++;
                                seq_dat[curRecSeq2].site[curSite] = '3';
                                seq_dat[curRecSeq2].baseCount3[curCodPos]++;
                                seq_dat[curRecSeq2].baseCount1[curCodPos]--;
                            }
                            else if (seq_dat[curRecSeq2].site[curSite] == '2')
                            {
                                /*site in first seq is 3 and in second is 2
                                 so swap and update counts*/
                                curCodPos=sitePosToCodPos(curSite);
                                seq_dat[curRecSeq1].site[curSite] ='2';
                                seq_dat[curRecSeq1].baseCount3[curCodPos]--;
                                seq_dat[curRecSeq1].baseCount2[curCodPos]++;
                                seq_dat[curRecSeq2].site[curSite] = '3';
                                seq_dat[curRecSeq2].baseCount3[curCodPos]++;
                                seq_dat[curRecSeq2].baseCount2[curCodPos]--;
                            }
                        }
                    }
                }
                /*skip one segment and go to the one after that*/
                curSeg=curSeg+2;
            }
        }
    }
    /* check basects after recombination */
    // 	for (i=1; i<=curPopSize; i++)
    // 	{
    // 		for (curCodPos = 1; curCodPos <= 3; curCodPos++)
    // 		{
    // 			baseCount0=0;
    // 			baseCount1=0;
    // 			for (j=curCodPos; j<=userPrefp->seqLen; j+=3)
    // 			{
    // 				if (seq_dat[i].site[j] == '0')
    // 					baseCount0++;
    // 				else
    // 					baseCount1++;
    // 			}
    // 			if (baseCount1 != seq_dat[i].baseCount1[curCodPos])
    // 			{
    // 				printf("basect mismatch after rec:\n");
    // 				printf("curSeq:                         %0ld\n", i);
    // 				printf("curCodPos:                      %0ld\n", curCodPos);
    // 				printf("baseCount1:                     %0ld\n", baseCount1);
    // 				printf("seq_dat[i].baseCount1[curCodPos]:   %0ld\n", seq_dat[i].baseCount1[curCodPos]);
    // 			}
    // 		}
    // 	}
    
    //int j;
    //for (curSeq=1; curSeq <= 1; curSeq++)
    //{
    //  for (j=1; j<=userPrefp->seqLen; j=j+3) {
    //    printf ("%c%c%c", seq_dat[curSeq].site[j], seq_dat[curSeq].site[j+1], seq_dat[curSeq].site[j+2]);
    //}
    //}
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*Function that evaluates the fitness of all individuals and then using multinomial
 sampling finds the numbers of offspring for each individual and creates the next
 generation of individuals. Function uses sArray to find fitness of individuals
 the following are initialized in initDat() */
//	double wData[MAXSEQNUM+1];				/*array that stores the weight for each sequence*/
//	double freqs[MAXSEQNUM+1];				/*expected frequencies of current individuals*/
//	long newFreqs[MAXSEQNUM+1];			/*obtained frequencies of offsprings for current generation*/
long getNextGenSeqs(struct userPref_str *userPrefp, long curGenNum)
{
    long curSeq;					/*refers to current sequence*/
    long curParent;				/*refers to the current parent sequence*/
    long curChild;					/*refers to all child sequences of a given parent*/
    long curSite;					/*refers to current site*/
    long newPopSize;				/*the new population size*/
    double wSum;							/*sum of wData*/
    void *tempPV;							/*temporary pointer variable for swapping sequence pointers*/
    
#ifdef DEBUGMODE
    printf("Finding next gen seqs\n");
    halt;
#endif
    
    
    
    //if (curGenNum <= userPrefp->uswitch) {
    //    population_size = userPrefp->initSeqNum;
    //}
    //else if (curGenNum > userPrefp->uswitch) {
    //    population_size = userPrefp->SeqNum_b;
    //}
    /*Initialize weights to 1.0 each*/
    for (curSeq=1; curSeq <= population_size; curSeq++)
    {
        wData[curSeq] = 1.0;
    }
    /*Loops through all variable sites and finds the weights for each individual*/
    for (curSite=1; curSite <= userPrefp->seqLen; curSite++)
    {
        if (polySite[curSite])
        {
            for (curSeq=1; curSeq<=population_size; curSeq++)
            {
                if (seq_dat[curSeq].site[curSite] != siteData[curSite].ancSite)
                {
                    wData[curSeq] *= (1.0 + sArray[curSite]);						/* multiplicative fitness implemented !!! */
                }
            }
        }
    }
    /*Finds the sum of weights to get relative weights*/
    wSum=0.0;
    for(curSeq=1;curSeq<=population_size;curSeq++)
    {
        wSum += wData[curSeq];
    }
    /*expected frequencies in the nextgeneration = relative weights*/
    for(curSeq=1;curSeq<=population_size;curSeq++)
    {
        freqs[curSeq]=wData[curSeq]/wSum;
    }
    /*the population size for a given generation is given in the popSizeArray*/
    newPopSize = popSizeArray[curGenNum+1];
    
    /*take multinomial deviates on the expected frequencies*/
    multdev(freqs, population_size, newPopSize, newFreqs);
    /*create new generation of individuals from the new frequences from
     multinomial sampling*/
    curSeq=1;
    for(curParent=1;curParent<=population_size;curParent++)
    {
        for(curChild=1;curChild<=newFreqs[curParent];curChild++)
        {
            //			new_seq_dat[curSeq] = seq_dat[curParent];			/* HAcheck */
            curSeq++;
        }
    }
    /*check if sampling is working correctly*/
    //if ((curSeq-1) != newPopSize)
    //	errorOut(("Invalid number of new sequences(%lu!=%lu) in getNextGenSeqs", curSeq, newPopSize));
    /*population size=population size in the new generation*/
    //	curPopSize=newPopSize;
    
    /*swap pointers*/
    tempPV = seq_dat;
    //	seq_dat = new_seq_dat;
    //	new_seq_dat = tempPV;
    
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*Function that evaluates the fitness of all individuals and then using multinomial
 sampling finds the numbers of offspring for each individual and creates the next
 generation of individuals. Function uses sArray to find fitness of individuals
 the following are initialized in initDat() */
//double wData[MAXSEQNUM+1];				/*array that stores the weight for each sequence*/
//double freqs[MAXSEQNUM+1];				/*expected frequencies of current individuals*/
//long newFreqs[MAXSEQNUM+1];			/*obtained frequencies of offsprings for current generation*/
long getNextGenSeqs_MCU(struct userPref_str *userPrefp, long curGenNum)
{
    int i, j;
    long curSeq, curCodPos;						/*refers to current sequence*/
    //	long curSite, curCodPos;
    long curParent;				/*refers to the current parent sequence*/
    long curChild;					/*refers to all child sequences of a given parent*/
    long selected_site_num;
    //	long curSite;					/*refers to current site*/
    long newPopSize;				/*the new population size*/
    double wSum;							/*sum of wData*/
    void *tempPV;							/*temporary pointer variable for swapping sequence pointers*/
    int cur_codon_index;
    int prev_codon_index;
    char change_type;
    double fitness;
    
#ifdef DEBUGMODE
    printf("Finding next gen seqs\n");
    halt;
#endif
    
    
    /* get fitness for each seq and calc ave */
    /* 150611matsumoto: fitness is given based on codon */
    
    /*for (i=0;i<=63;i++) {
     if (check==0) {
     cod_w_list[i] = cod_w_list_A[i];
     }
     else if (check==1 && (curGenNum <= userPrefp->uswitch)) {
     cod_w_list[i] = cod_w_list_B[i];
     }
     else if (check==1 && (curGenNum > userPrefp->uswitch)) {
     cod_w_list[i] = cod_w_list_C[i];
     }
     }*/
    
    wSum=0.0;
    
    for (curSeq=1; curSeq <= population_size; curSeq++)
    {
        wData[curSeq]=1.00;
        for (j=1; j<=userPrefp->seqLen; j=j+3) {
            if (seq_dat[curSeq].site[j+1] == '0') {
                if (seq_dat[curSeq].site[j] == '0') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=0;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=1;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=2;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=3;
                    }
                }
                else if (seq_dat[curSeq].site[j] == '1') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=4;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=5;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=6;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=7;
                    }
                }
                else if (seq_dat[curSeq].site[j] == '2') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=8;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=9;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=10;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=11;
                    }
                }
                else if (seq_dat[curSeq].site[j] == '3') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=12;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=13;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=14;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=15;
                    }
                }
            }
            
            if (seq_dat[curSeq].site[j+1] == '1') {
                if (seq_dat[curSeq].site[j] == '0') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=16;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=17;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=18;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=19;
                    }
                }
                else if (seq_dat[curSeq].site[j] == '1') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=20;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=21;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=22;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=23;
                    }
                }
                else if (seq_dat[curSeq].site[j] == '2') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=24;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=25;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=26;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=27;
                    }
                }
                else if (seq_dat[curSeq].site[j] == '3') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=28;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=29;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=30;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=31;
                    }
                }
            }
            
            if (seq_dat[curSeq].site[j+1] == '2') {
                if (seq_dat[curSeq].site[j] == '0') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=32;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=33;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=34;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=35;
                    }
                }
                else if (seq_dat[curSeq].site[j] == '1') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=36;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=37;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=38;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=39;
                    }
                }
                else if (seq_dat[curSeq].site[j] == '2') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=40;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=41;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=42;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=43;
                    }
                }
                else if (seq_dat[curSeq].site[j] == '3') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=44;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=45;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=46;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=47;
                    }
                }
            }
            
            if (seq_dat[curSeq].site[j+1] == '3') {
                if (seq_dat[curSeq].site[j] == '0') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=48;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=49;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=50;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=51;
                    }
                }
                else if (seq_dat[curSeq].site[j] == '1') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=52;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=53;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=54;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=55;
                    }
                }
                else if (seq_dat[curSeq].site[j] == '2') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=56;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=57;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=58;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=59;
                    }
                }
                else if (seq_dat[curSeq].site[j] == '3') {
                    if (seq_dat[curSeq].site[j+2] == '0') {
                        cur_codon_index=60;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '1') {
                        cur_codon_index=61;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '2') {
                        cur_codon_index=62;
                    }
                    else if (seq_dat[curSeq].site[j+2] == '3') {
                        cur_codon_index=63;
                    }
                }
            }
            
            if (prev_seq_dat[1].site[j+1] == '0') {
                if (prev_seq_dat[1].site[j] == '0') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=0;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=1;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=2;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=3;
                    }
                }
                else if (prev_seq_dat[1].site[j] == '1') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=4;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=5;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=6;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=7;
                    }
                }
                else if (prev_seq_dat[1].site[j] == '2') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=8;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=9;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=10;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=11;
                    }
                }
                else if (prev_seq_dat[1].site[j] == '3') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=12;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=13;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=14;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=15;
                    }
                }
            }
            
            if (prev_seq_dat[1].site[j+1] == '1') {
                if (prev_seq_dat[1].site[j] == '0') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=16;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=17;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=18;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=19;
                    }
                }
                else if (prev_seq_dat[1].site[j] == '1') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=20;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=21;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=22;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=23;
                    }
                }
                else if (prev_seq_dat[1].site[j] == '2') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=24;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=25;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=26;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=27;
                    }
                }
                else if (prev_seq_dat[1].site[j] == '3') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=28;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=29;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=30;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=31;
                    }
                }
            }
            
            if (prev_seq_dat[1].site[j+1] == '2') {
                if (prev_seq_dat[1].site[j] == '0') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=32;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=33;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=34;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=35;
                    }
                }
                else if (prev_seq_dat[1].site[j] == '1') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=36;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=37;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=38;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=39;
                    }
                }
                else if (prev_seq_dat[1].site[j] == '2') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=40;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=41;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=42;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=43;
                    }
                }
                else if (prev_seq_dat[1].site[j] == '3') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=44;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=45;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=46;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=47;
                    }
                }
            }
            
            if (prev_seq_dat[1].site[j+1] == '3') {
                if (prev_seq_dat[1].site[j] == '0') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=48;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=49;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=50;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=51;
                    }
                }
                else if (prev_seq_dat[1].site[j] == '1') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=52;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=53;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=54;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=55;
                    }
                }
                else if (prev_seq_dat[1].site[j] == '2') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=56;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=57;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=58;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=59;
                    }
                }
                else if (prev_seq_dat[1].site[j] == '3') {
                    if (prev_seq_dat[1].site[j+2] == '0') {
                        prev_codon_index=60;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '1') {
                        prev_codon_index=61;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '2') {
                        prev_codon_index=62;
                    }
                    else if (prev_seq_dat[1].site[j+2] == '3') {
                        prev_codon_index=63;
                    }
                }
            }
            
            
            /* 151221 matsumoto: for synonymous site sequence, "3", "2", "1", "0" depending  on the number of prefered codon */
            /* 150617 matsumoto: change_type is "s", "r", ""p", "u", "n": stop is stop codon, r is replacement, u is silent (p->u 3rd pos), p is silent (u->p 3rd pos), n is silent (p->p or u->u 3rd pos) */
            change_type = userPrefp->codon_change[prev_codon_index][cur_codon_index];
            //if (prev_codon_index != cur_codon_index) {
            //    printf ("%d, %d, %c\n", prev_codon_index, cur_codon_index, userPrefp->codon_change[prev_codon_index][cur_codon_index]);
            //}
            
            if (change_type == '3') {
                fitness = (userPrefp->fitness_p) * (userPrefp->fitness_p) * (userPrefp->fitness_p);
            }
            else if (change_type == '2') {
                fitness = (userPrefp->fitness_u) * (userPrefp->fitness_p) * (userPrefp->fitness_p);
            }
            else if (change_type == '1') {
                fitness = (userPrefp->fitness_u) * (userPrefp->fitness_u) * (userPrefp->fitness_p);
            }
            else if (change_type == '0') {
                fitness = (userPrefp->fitness_u) * (userPrefp->fitness_u) * (userPrefp->fitness_u);
            }
            
            wData[curSeq] = wData[curSeq] * fitness;
        }
        
        
        wSum += wData[curSeq];
    }
    
    /*expected frequencies in the nextgeneration = relative weights*/
    for (curSeq=1; curSeq<=population_size; curSeq++)
    {
        //printf ("%f\n", wData[curSeq]);
        freqs[curSeq] = wData[curSeq] / wSum;
    }
    
    /*the population size for a given generation is given in the popSizeArray*/
    newPopSize = popSizeArray[curGenNum+1];
    //printf ("%ld\t%ld\n", newPopSize, population_size);
    
    
    /*take multinomial deviates on the expected frequencies*/
    multdev(freqs, population_size, newPopSize, newFreqs);
    
    curSeq=0;
    for (curParent=1; curParent <= population_size; curParent++)
    {
        for (curChild=1; curChild <= newFreqs[curParent]; curChild++)
        {
            curSeq++;
        }
    }
    //printf ("%d\t%d\n", curParent, curSeq);
    /*check if sampling is working correctly*/
    //if (curSeq != population_size)
    //	errorOut(("Invalid number of new sequences(%lu!=%lu) in getNextGenSeqs_MCU", curSeq, population_size));
    
    /* option 1 - anoop's original pointer method */
    /*create new generation of individuals from the new frequences from multinomial sampling*/
    curSeq=0;
    
    /* 150602 matsumoto: population pass from prerun to initrun */
    if (check==0) {
        if (curGenNum >= userPrefp->initRunGen) {
            population_size = userPrefp->initSeqNum;
        }
    }
    for (curParent=1; curParent<=population_size; curParent++)
    {
        for (curChild=1; curChild <= newFreqs[curParent]; curChild++)
        {
            curSeq++;
            temp_seq_dat[curSeq] = seq_dat[curParent];
        }
    }
    
    /*swap pointers*/
    tempPV = seq_dat;
    seq_dat = temp_seq_dat;
    temp_seq_dat = tempPV;
    
    /* end option 1 */
    
    
    /* HA pointer option 2 - this is very slow */
    /* copy seqdat to temp_seq_dat */
    // 	for (curSeq=1; curSeq<= userPrefp->initSeqNum; curSeq++)
    // 	{
    //  		for (curSite = 1; curSite <= userPrefp->seqLen; curSite++)
    //  			temp_seq_dat[curSeq].site[curSite] = seq_dat[curSeq].site[curSite];
    // 		for (curCodPos = 1; curCodPos <= 3; curCodPos++)
    // 		{
    // 			 temp_seq_dat[curSeq].baseCount0[curCodPos] = seq_dat[curSeq].baseCount0[curCodPos];
    // 			 temp_seq_dat[curSeq].baseCount1[curCodPos] = seq_dat[curSeq].baseCount1[curCodPos];
    // 		}
    //  		temp_seq_dat[curSeq].tot_1_ct = seq_dat[curSeq].tot_1_ct;
    // 	}
    
    /*create new generation of individuals from the new frequences from multinomial sampling*/
    // 	curSeq=0;
    // 	for (curParent=1; curParent <= curPopSize; curParent++)
    // 	{
    // 		for (curChild=1; curChild <= newFreqs[curParent]; curChild++)
    // 		{
    // 			curSeq++;
    // 			/* copy seqdat from temp_seq_dat */
    // 			for (curSeq=1; curSeq<= userPrefp->initSeqNum; curSeq++)
    // 			{
    // 				for (curSite = 1; curSite <= userPrefp->seqLen; curSite++)
    // 					seq_dat[curSeq].site[curSite] = temp_seq_dat[curParent].site[curSite];
    // 				for (curCodPos = 1; curCodPos <= 3; curCodPos++)
    // 				{
    // 					 seq_dat[curSeq].baseCount0[curCodPos] = temp_seq_dat[curParent].baseCount0[curCodPos];
    // 					 seq_dat[curSeq].baseCount1[curCodPos] = temp_seq_dat[curParent].baseCount1[curCodPos];
    // 				}
    // 				seq_dat[curSeq].tot_1_ct = temp_seq_dat[curSeq].tot_1_ct;
    // 			}
    // 		}
    // 	}
    
    /*population size=population size in the new generation*/
    //	curPopSize=newPopSize;
    
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*During reproduction sequences may get lost due to multinomial sampling and/or selection.
 Polymorphic sites may become monomorphic or the derived state might get fixed
 in the generation. This function updates the siteData structure and the polySites
 flag. Also stores the fixation positions, svals and ancsite so that these can be
 stored in the outputfile
 */
long updateCounts(struct userPref_str *userPrefp, long curGenNum)
{
    long curSite;					/*refers to the current site*/
    long curSeq;					/*refers to the current sequence*/
    long baseCount0;				/*number of 0 bases at curSite across all seqs*/
    long baseCount1;				/*number of 1 bases at curSite across all seqs*/
    long baseCount2;				/*number of 2 bases at curSite across all seqs*/
    long baseCount3;				/*number of 3 bases at curSite across all seqs*/
    
#ifdef DEBUGMODE
    printf("Updating Counts\n");
    halt;
#endif
    
    //if (curGenNum <= userPrefp->uswitch) {
    //    population_size = userPrefp->initSeqNum;
    //}
    //else if (curGenNum > userPrefp->uswitch) {
    //    population_size = userPrefp->SeqNum_b;
    //}
    /*Loop across all the sites*/
    for(curSite=1; curSite <= userPrefp->seqLen; curSite++)
    {
        /*check only if the polySite flag is set*/
        if(polySite[curSite] == 1)
        {
            /*Find the frequencies of 0s and 1s at the site*/
            baseCount0=0;
            baseCount1=0;
            baseCount2=0;
            baseCount3=0;
            
            for (curSeq=1; curSeq<=population_size; curSeq++)
            {
                if (seq_dat[curSeq].site[curSite] == '0') {
                    baseCount0++;
                }
                else if (seq_dat[curSeq].site[curSite] == '1') {
                    baseCount1++;
                }
                else if (seq_dat[curSeq].site[curSite] == '2') {
                    baseCount2++;
                }
                else if (seq_dat[curSeq].site[curSite] == '3') {
                    baseCount3++;
                }
            }
            /*update the counts in the siteData structure*/
            siteData[curSite].freqs0 = baseCount0;
            siteData[curSite].freqs1 = baseCount1;
            siteData[curSite].freqs2 = baseCount2;
            siteData[curSite].freqs3 = baseCount3;
            
            /*Check if any of the basecounts is 0*/
            if((baseCount1==0) && (baseCount2==0) &&(baseCount3==0))
            {
                /*if so check if the ancestral site is same as the base whose count is zero
                 if so a fixation. here check if 0 has been fixed*/
                if(siteData[curSite].ancSite == '1')
                {
                    /*if fixation update counts in the siteData structure*/
                    siteData[curSite].ancSite = '0';
                    if (check==1) {
                        siteData[curSite].fixNum10++;
                        /*store the svalues positions and fixed site to be output in the writerepdata*/
                        popFixNum++;
                        fixSite[popFixNum] = '0';
                        fixSitePos[popFixNum] = curSite;
                        fixSiteSval[popFixNum] = sArray[curSite];
                        /*at every fixation resample for the sArray*/
                        sArray[curSite] = 0.0;
                    }
                }
                else if(siteData[curSite].ancSite == '2')
                {
                    /*if fixation update counts in the siteData structure*/
                    siteData[curSite].ancSite = '0';
                    if (check==1) {
                        siteData[curSite].fixNum20++;
                        /*store the svalues positions and fixed site to be output in the writerepdata*/
                        popFixNum++;
                        fixSite[popFixNum] = '0';
                        fixSitePos[popFixNum] = curSite;
                        fixSiteSval[popFixNum] = sArray[curSite];
                        /*at every fixation resample for the sArray*/
                        sArray[curSite] = 0.0;
                    }
                }
                else if(siteData[curSite].ancSite == '3')
                {
                    /*if fixation update counts in the siteData structure*/
                    siteData[curSite].ancSite = '0';
                    if (check==1) {
                        siteData[curSite].fixNum30++;
                        /*store the svalues positions and fixed site to be output in the writerepdata*/
                        popFixNum++;
                        fixSite[popFixNum] = '0';
                        fixSitePos[popFixNum] = curSite;
                        fixSiteSval[popFixNum] = sArray[curSite];
                        /*at every fixation resample for the sArray*/
                        sArray[curSite] = 0.0;
                    }
                }
                polySite[curSite] = 0;
            }
            
            else if((baseCount0==0) && (baseCount2==0) &&(baseCount3==0))
            {
                /*if so check if the ancestral site is same as the base whose count is zero
                 if so a fixation. here check if 0 has been fixed*/
                if(siteData[curSite].ancSite == '0')
                {
                    /*if fixation update counts in the siteData structure*/
                    siteData[curSite].ancSite = '1';
                    if (check==1) {
                        siteData[curSite].fixNum01++;
                        /*store the svalues positions and fixed site to be output in the writerepdata*/
                        popFixNum++;
                        fixSite[popFixNum] = '1';
                        fixSitePos[popFixNum] = curSite;
                        fixSiteSval[popFixNum] = sArray[curSite];
                        /*at every fixation resample for the sArray*/
                        sArray[curSite] = 0.0;
                    }
                }
                else if(siteData[curSite].ancSite == '2')
                {
                    /*if fixation update counts in the siteData structure*/
                    siteData[curSite].ancSite = '1';
                    if (check==1) {
                        siteData[curSite].fixNum21++;
                        /*store the svalues positions and fixed site to be output in the writerepdata*/
                        popFixNum++;
                        fixSite[popFixNum] = '1';
                        fixSitePos[popFixNum] = curSite;
                        fixSiteSval[popFixNum] = sArray[curSite];
                        /*at every fixation resample for the sArray*/
                        sArray[curSite] = 0.0;
                    }
                }
                else if(siteData[curSite].ancSite == '3')
                {
                    /*if fixation update counts in the siteData structure*/
                    siteData[curSite].ancSite = '1';
                    if (check==1) {
                        siteData[curSite].fixNum31++;
                        /*store the svalues positions and fixed site to be output in the writerepdata*/
                        popFixNum++;
                        fixSite[popFixNum] = '1';
                        fixSitePos[popFixNum] = curSite;
                        fixSiteSval[popFixNum] = sArray[curSite];
                        /*at every fixation resample for the sArray*/
                        sArray[curSite] = 0.0;
                    }
                }
                polySite[curSite] = 0;
            }
            
            else if((baseCount0==0) && (baseCount1==0) &&(baseCount3==0))
            {
                /*if so check if the ancestral site is same as the base whose count is zero
                 if so a fixation. here check if 0 has been fixed*/
                if(siteData[curSite].ancSite == '0')
                {
                    /*if fixation update counts in the siteData structure*/
                    siteData[curSite].ancSite = '2';
                    if (check==1) {
                        siteData[curSite].fixNum02++;
                        /*store the svalues positions and fixed site to be output in the writerepdata*/
                        popFixNum++;
                        fixSite[popFixNum] = '2';
                        fixSitePos[popFixNum] = curSite;
                        fixSiteSval[popFixNum] = sArray[curSite];
                        /*at every fixation resample for the sArray*/
                        sArray[curSite] = 0.0;
                    }
                }
                else if(siteData[curSite].ancSite == '1')
                {
                    /*if fixation update counts in the siteData structure*/
                    siteData[curSite].ancSite = '2';
                    if (check==1) {
                        siteData[curSite].fixNum12++;
                        /*store the svalues positions and fixed site to be output in the writerepdata*/
                        popFixNum++;
                        fixSite[popFixNum] = '2';
                        fixSitePos[popFixNum] = curSite;
                        fixSiteSval[popFixNum] = sArray[curSite];
                        /*at every fixation resample for the sArray*/
                        sArray[curSite] = 0.0;
                    }
                }
                else if(siteData[curSite].ancSite == '3')
                {
                    /*if fixation update counts in the siteData structure*/
                    siteData[curSite].ancSite = '2';
                    if (check==1) {
                        siteData[curSite].fixNum32++;
                        /*store the svalues positions and fixed site to be output in the writerepdata*/
                        popFixNum++;
                        fixSite[popFixNum] = '2';
                        fixSitePos[popFixNum] = curSite;
                        fixSiteSval[popFixNum] = sArray[curSite];
                        /*at every fixation resample for the sArray*/
                        sArray[curSite] = 0.0;
                    }
                }
                polySite[curSite] = 0;
            }
            
            else if((baseCount0==0) && (baseCount1==0) &&(baseCount2==0))
            {
                /*if so check if the ancestral site is same as the base whose count is zero
                 if so a fixation. here check if 0 has been fixed*/
                if(siteData[curSite].ancSite == '0')
                {
                    /*if fixation update counts in the siteData structure*/
                    siteData[curSite].ancSite = '3';
                    if (check==1) {
                        siteData[curSite].fixNum03++;
                        /*store the svalues positions and fixed site to be output in the writerepdata*/
                        popFixNum++;
                        fixSite[popFixNum] = '3';
                        fixSitePos[popFixNum] = curSite;
                        fixSiteSval[popFixNum] = sArray[curSite];
                        /*at every fixation resample for the sArray*/
                        sArray[curSite] = 0.0;
                    }
                }
                else if(siteData[curSite].ancSite == '1')
                {
                    /*if fixation update counts in the siteData structure*/
                    siteData[curSite].ancSite = '3';
                    if (check==1) {
                        siteData[curSite].fixNum13++;
                        /*store the svalues positions and fixed site to be output in the writerepdata*/
                        popFixNum++;
                        fixSite[popFixNum] = '3';
                        fixSitePos[popFixNum] = curSite;
                        fixSiteSval[popFixNum] = sArray[curSite];
                        /*at every fixation resample for the sArray*/
                        sArray[curSite] = 0.0;
                    }
                }
                else if(siteData[curSite].ancSite == '2')
                {
                    /*if fixation update counts in the siteData structure*/
                    siteData[curSite].ancSite = '3';
                    if (check==1) {
                        siteData[curSite].fixNum23++;
                        /*store the svalues positions and fixed site to be output in the writerepdata*/
                        popFixNum++;
                        fixSite[popFixNum] = '3';
                        fixSitePos[popFixNum] = curSite;
                        fixSiteSval[popFixNum] = sArray[curSite];
                        /*at every fixation resample for the sArray*/
                        sArray[curSite] = 0.0;
                    }
                }
                polySite[curSite] = 0;
            }
        }
    }
    /*If number of fixations greater than maximum permitted then error*/
    if(popFixNum>MAXFIXATIONS)
    errorOut(("Exceeded maximum fixations permitted, redefine MAXFIXATIONS"));
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*
 Function used to fill and modify the sArray. This function is called
 when initializing and when fixations occur
 
 */
double sampleSval(struct userPref_str *userPrefp, long sitePos, char ancSite)
{
    long codPos;
    //	const double Nes = 10.0;
    const double Nes = 0.00001;
    const double Ne = 500.0;
    const double s01 = Nes/Ne;
    const double s02 = Nes/Ne;
    const double s03 = Nes/Ne;
    const double s11 = Nes/Ne;
    const double s12 = Nes/Ne;
    const double s13 = Nes/Ne;
    double prob = 1.0/10.0;					/*the fraction of times the function returns a nonzero value*/
    /*changed the variables to consts - AJ 07/21/03*/
    /* HA changed to neutral for now 121111 */
    
    return 0.0;
    
    
    prob = 1.02;		/* all mutations will be neutral */
    if(ran2() > prob)
    return 0.0;
    codPos=sitePosToCodPos(sitePos);
    switch(codPos)
    {
        case	1:
            if(ancSite == '1')
            return s11;
            else
            return s01;
            break;
        case	2:
            if(ancSite == '1')
            return s12;
            else
            return s02;
            break;
        case	3:
            if(ancSite == '1')
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
 information, fixation information including position, fixed site and Svals, segregating
 sites information and the common ancestor of the samples are also printed ou*/
long writeRepData(struct userPref_str *userPrefp, long curGenNum)
{
    long i, j, k, sum_S, sum_R, sum, sum01, sum02, sum03, sum10, sum12, sum13, sum20, sum21, sum23, sum30, sum31, sum32;						/*for loop*/
    long numPopPeriods;			/*number of periods for the population fluctuation scenario*/
    long startPeriod;				/*starting generation of a population period*/
    long endPeriod;				/*ending generation of a population period*/
    long curSeq;					/*refers to the current sampled sequence*/
    double ave_MCU_silent;
    double ave_MCU_replacement;
    double ave_MCU;
    
    /*sample sequences and collect data from the samples*/
    //getRepData(userPrefp);
    
    //if (curGenNum <= userPrefp->uswitch) {
    //    population_size = userPrefp->initSeqNum;
    //}
    //else if (curGenNum > userPrefp->uswitch) {
    //    population_size = userPrefp->SeqNum_b;
    //}
    
    //	printf("\tWriting Replicate information to file\n");
    /*Outputs the counts obtained during the current replicate*/
    fprintf(fpOut, "Replicate:%6ld\n", curRepNum);
    fprintf(fpOut, "--- seq_dat[curSeq].tot_1_ct -------------\n");
    //	sum_S = 0;
    //    sum_R = 0;
    //	for(curSeq=1; curSeq<= population_size; curSeq++)
    //	{
    //		fprintf(fpOut, "%0ld\n", seq_dat[curSeq].tot_1_ct);
    //		sum_S += seq_dat[curSeq].tot_1_ct_silent;
    //        sum_R += seq_dat[curSeq].tot_1_ct_replacement;
    //	}
    //	ave_MCU_silent = (double) sum_S / (double) population_size;
    //    ave_MCU_replacement = (double) sum_R / (double) population_size;
    //	fprintf(fpOut, "ave_silent: %0.4f\n", ave_MCU_silent);
    //    fprintf(fpOut, "ave_replacement: %0.4f\n", ave_MCU_replacement);
    fprintf(fpOut, "----------------\n");
    
    
    /*prints out the sitedata information plus sArray*/
    fprintf(fpOut, "SiteData(ancSite:sVal:freqs0:freqs1:freqs2:freqs3:fix01:fix02:fix03:fix10:fix12:fix13:fix20:fix21:fix23:fix30:fix31:fix32)\n");
    fprintf(fpOut, "----------------------------------------------------------------\n");
    for(i=1; i<= userPrefp->seqLen; i++)
    {
        fprintf(fpOut, "%c:", siteData[i].ancSite);
        fprintf(fpOut, "%.10f:", sArray[i]);
        fprintf(fpOut, "%lu:", siteData[i].freqs0);
        fprintf(fpOut, "%lu:", siteData[i].freqs1);
        fprintf(fpOut, "%lu:", siteData[i].freqs2);
        fprintf(fpOut, "%lu:", siteData[i].freqs3);
        fprintf(fpOut, "%lu:", siteData[i].fixNum01);
        fprintf(fpOut, "%lu:", siteData[i].fixNum02);
        fprintf(fpOut, "%lu:", siteData[i].fixNum03);
        fprintf(fpOut, "%lu:", siteData[i].fixNum10);
        fprintf(fpOut, "%lu:", siteData[i].fixNum12);
        fprintf(fpOut, "%lu:", siteData[i].fixNum13);
        fprintf(fpOut, "%lu:", siteData[i].fixNum20);
        fprintf(fpOut, "%lu:", siteData[i].fixNum21);
        fprintf(fpOut, "%lu:", siteData[i].fixNum23);
        fprintf(fpOut, "%lu:", siteData[i].fixNum30);
        fprintf(fpOut, "%lu:", siteData[i].fixNum31);
        fprintf(fpOut, "%lu:", siteData[i].fixNum32);
        
    }
    fprintf(fpOut, "\n");
    fprintf(fpOut, "\n");
    
    
    /* prints out the sequences of all individuals */
    fprintf(fpOut, "whole sequence of each sample\n");
    for(i=1;i<=userPrefp->SeqNum_b;i++)
    {
        fprintf(fpOut, "seq_%ld_							= \"", i);
        for(j=1;j<=userPrefp->seqLen;j++)
        {
            fprintf(fpOut, "%c", seq_dat[i].site[j]);
        }
        fprintf(fpOut, "\"\n");
    }
    
    fprintf(fpOut, "ancestor_							= \"");
    for(j=1;j<=userPrefp->seqLen;j++)
    {
        fprintf(fpOut, "%c", siteData[j].ancSite);
    }
    fprintf(fpOut, "\"\n");
    
    /*fprintf(fpOut, "\n");
     fprintf(fpOut, "segregating replacement sites\n");
     for(i=1;i<=userPrefp->sampleNum;i++)
     {
     curSeq=sampleSeqs[i];
     for(j=1;j<=sampleSegNum;j++)
     {
     if (sampleSegPos[j] % 2 != 0) {
     fprintf(fpOut, "%c", seq_dat[curSeq].site[sampleSegPos[j]]);
     }
     }
     fprintf(fpOut, "\n");
     }
     fprintf(fpOut, "\n");
     fprintf(fpOut, "invariant sites\n");
     fprintf(fpOut, "------------\n");
     
     j=1;
     for (k=1;k<=userPrefp->seqLen;k++) {
     if (k!=sampleSegPos[j]) {
     fprintf(fpOut, "%lu ", k);
     }
     if (k==sampleSegPos[j]) {
     j++;
     }
     }
     
     fprintf(fpOut, "\n");
     fprintf(fpOut, "\n");
     fprintf(fpOut, "invariant silent sites\n");
     j=1;
     for (k=1;k<=userPrefp->seqLen;k++) {
     if (k!=sampleSegPos[j] && k % 2 == 0) {
     fprintf(fpOut, "%c", seq_dat[1].site[k]);
     }
     if (k==sampleSegPos[j]) {
     j++;
     }
     }
     fprintf(fpOut, "\n");
     fprintf(fpOut, "\n");
     fprintf(fpOut, "invariant replacement sites\n");
     j=1;
     for (k=1;k<=userPrefp->seqLen;k++) {
     if (k!=sampleSegPos[j] && k % 2 != 0) {
     fprintf(fpOut, "%c", seq_dat[1].site[k]);
     }
     if (k==sampleSegPos[j]) {
     j++;
     }
     }*/
    
    //fprintf(fpOut, "//\n");
    //fflush(fpOut);
    //fflush(fpOut2);
    //fflush(fpOut_segsites);
    //fflush(fpOut_MCU);
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*Collects the information from each replicate. This involves picking a sample set
 from the population and collecting information from the sample.*/
long getRepData(struct userPref_str *userPrefp)
{
    long curSite;					/*refers to the current site*/
    long freq0;					/*frequencies of 0 at the current site*/
    long freq1;					/*frequencies of 1 at the current site*/
    long freq2;					/*frequencies of 2 at the current site*/
    long freq3;					/*frequencies of 3 at the current site*/
    long curSample;				/*refers to the current sample num*/
    
    //	printf("\tGetting Replicate Information\n");
    /*reset counts*/
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
    
    /*picks sampleNum sequences*/
    getSample(userPrefp, curGenNum);
    for (curSite=1; curSite<= userPrefp->seqLen; curSite++)
    {
        /*if the site is not polymorphic in the population then no need to process it*/
        if (polySite[curSite] == 1)
        {
            freq0=freq1=freq2=freq3=0;
            for(curSample=1; curSample <= userPrefp->sampleNum; curSample++)
            {
                if(seq_dat[sampleSeqs[curSample]].site[curSite] == '0') {
                    freq0++;
                }
                else if(seq_dat[sampleSeqs[curSample]].site[curSite] == '1') {
                    freq1++;
                }
                else if(seq_dat[sampleSeqs[curSample]].site[curSite] == '2') {
                    freq2++;
                }
                else if(seq_dat[sampleSeqs[curSample]].site[curSite] == '3') {
                    freq3++;
                }
            }
            /*if site is polymorphic in sample*/
            if ((freq1!=userPrefp->sampleNum)&&(freq0!=userPrefp->sampleNum)&&(freq2!=userPrefp->sampleNum)&&(freq3!=userPrefp->sampleNum))
            {
                /*depending on the ancestral state update the counters*/
                if (siteData[curSite].ancSite == '1')
                {
                    /*150605 matsumoto: we do not allow multiple mutation at the same site, so there are two states for each polymorphic site*/
                    if (freq0 > 0) {
                        sampleSeg10[sitePosToCodPos(curSite)]++;
                    }
                    else if (freq2 > 0) {
                        sampleSeg12[sitePosToCodPos(curSite)]++;
                    }
                    else if (freq3 > 0) {
                        sampleSeg13[sitePosToCodPos(curSite)]++;
                    }
                }
                if (siteData[curSite].ancSite == '0')
                {
                    /*150605 matsumoto: we do not allow multiple mutation at the same site, so there are two states for each polymorphic site*/
                    if (freq1 > 0) {
                        sampleSeg01[sitePosToCodPos(curSite)]++;
                    }
                    else if (freq2 > 0) {
                        sampleSeg02[sitePosToCodPos(curSite)]++;
                    }
                    else if (freq3 > 0) {
                        sampleSeg03[sitePosToCodPos(curSite)]++;
                    }
                }
                if (siteData[curSite].ancSite == '2')
                {
                    /*150605 matsumoto: we do not allow multiple mutation at the same site, so there are two states for each polymorphic site*/
                    if (freq0 > 0) {
                        sampleSeg20[sitePosToCodPos(curSite)]++;
                    }
                    else if (freq1 > 0) {
                        sampleSeg21[sitePosToCodPos(curSite)]++;
                    }
                    else if (freq3 > 0) {
                        sampleSeg23[sitePosToCodPos(curSite)]++;
                    }
                }
                if (siteData[curSite].ancSite == '3')
                {
                    /*150605 matsumoto: we do not allow multiple mutation at the same site, so there are two states for each polymorphic site*/
                    if (freq0 > 0) {
                        sampleSeg30[sitePosToCodPos(curSite)]++;
                    }
                    else if (freq1 > 0) {
                        sampleSeg31[sitePosToCodPos(curSite)]++;
                    }
                    else if (freq2 > 0) {
                        sampleSeg32[sitePosToCodPos(curSite)]++;
                    }
                }
                /*keep track of the segregating sites*/
                sampleSegNum++;
                sampleSegPos[sampleSegNum]=curSite;
            }
            else
            {
                /*fixation of 0 due to sampling*/
                if (freq0==userPrefp->sampleNum) {
                    if (siteData[curSite].ancSite == '1') {
                        sampleFix10[sitePosToCodPos(curSite)]++;
                        sampleFixNum++;
                        sampleFixPos[sampleFixNum]=curSite;
                    }
                    else if (siteData[curSite].ancSite == '2') {
                        sampleFix20[sitePosToCodPos(curSite)]++;
                        sampleFixNum++;
                        sampleFixPos[sampleFixNum]=curSite;
                    }
                    else if (siteData[curSite].ancSite == '3') {
                        sampleFix30[sitePosToCodPos(curSite)]++;
                        sampleFixNum++;
                        sampleFixPos[sampleFixNum]=curSite;
                    }
                }
                
                else if (freq1==userPrefp->sampleNum) {
                    if (siteData[curSite].ancSite == '0') {
                        sampleFix01[sitePosToCodPos(curSite)]++;
                        sampleFixNum++;
                        sampleFixPos[sampleFixNum]=curSite;
                    }
                    else if (siteData[curSite].ancSite == '2') {
                        sampleFix21[sitePosToCodPos(curSite)]++;
                        sampleFixNum++;
                        sampleFixPos[sampleFixNum]=curSite;
                    }
                    else if (siteData[curSite].ancSite == '3') {
                        sampleFix31[sitePosToCodPos(curSite)]++;
                        sampleFixNum++;
                        sampleFixPos[sampleFixNum]=curSite;
                    }
                }
                
                else if (freq2==userPrefp->sampleNum) {
                    if (siteData[curSite].ancSite == '0') {
                        sampleFix02[sitePosToCodPos(curSite)]++;
                        sampleFixNum++;
                        sampleFixPos[sampleFixNum]=curSite;
                    }
                    else if (siteData[curSite].ancSite == '1') {
                        sampleFix12[sitePosToCodPos(curSite)]++;
                        sampleFixNum++;
                        sampleFixPos[sampleFixNum]=curSite;
                    }
                    else if (siteData[curSite].ancSite == '3') {
                        sampleFix32[sitePosToCodPos(curSite)]++;
                        sampleFixNum++;
                        sampleFixPos[sampleFixNum]=curSite;
                    }
                }
                
                else if (freq3==userPrefp->sampleNum) {
                    if (siteData[curSite].ancSite == '0') {
                        sampleFix03[sitePosToCodPos(curSite)]++;
                        sampleFixNum++;
                        sampleFixPos[sampleFixNum]=curSite;
                    }
                    else if (siteData[curSite].ancSite == '1') {
                        sampleFix13[sitePosToCodPos(curSite)]++;
                        sampleFixNum++;
                        sampleFixPos[sampleFixNum]=curSite;
                    }
                    else if (siteData[curSite].ancSite == '2') {
                        sampleFix23[sitePosToCodPos(curSite)]++;
                        sampleFixNum++;
                        sampleFixPos[sampleFixNum]=curSite;
                    }
                }
            }
        }
    }
    repTime=getTime(4);
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*randomly picks the current samplenum sequences from the population*/
long getSample(struct userPref_str *userPrefp, long curGenNum)
{
    long curSample;					/*refers to the current sample*/
    long randSeq;						/*randomly selected sequence*/
    long curSeq;
    long i;						/*loop through all samples*/
    
    //if (curGenNum <= userPrefp->uswitch) {
    //    population_size = userPrefp->initSeqNum;
    //}
    //else if (curGenNum > userPrefp->uswitch) {
    //    population_size = userPrefp->SeqNum_b;
    //}
    /*repeat loop until all random sequences are found*/
    curSample=0;
    while(curSample < userPrefp->sampleNum)
    {
        /*pick a random sequence*/
        randSeq=getRandLong_(1, population_size);
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
            for (i=0;i<=userPrefp->seqLen-1;i++) {
                printf(fpOut, "%d",seq_dat[randSeq]);
            }
        }
    }
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*function for freeing the memory allocated and closing all files*/
long endSim(struct userPref_str *userPrefp)
{
    free(seq_dat);
    free(temp_seq_dat);
    fclose(fpOut);
    fclose(fpOut2);
    fclose(fpOut_segsites);
    fclose(fpOut_MCU);
    return 0;
}
/*--------------------------------------------------------------------------------*/
/*function to check basecounts*/
long checkBaseCounts(struct userPref_str *userPrefp, long pos, long curGenNum)
{
    long i;
    long curSeq;
    long curSite;
    long baseCount0[3+1];
    long baseCount1[3+1];
    long baseCount2[3+1];
    long baseCount3[3+1];
    
    //if (curGenNum <= userPrefp->uswitch) {
    //    population_size = userPrefp->initSeqNum;
    //}
    //else if (curGenNum > userPrefp->uswitch) {
    //    population_size = userPrefp->SeqNum_b;
    //}
    
    for (curSeq=1; curSeq<=population_size; curSeq++)
    {
        baseCount0[1]=baseCount0[2]=baseCount0[3]=0;
        baseCount1[1]=baseCount1[2]=baseCount1[3]=0;
        baseCount2[1]=baseCount2[2]=baseCount2[3]=0;
        baseCount3[1]=baseCount3[2]=baseCount3[3]=0;
        
        for(curSite=1; curSite <= userPrefp->seqLen; curSite++)
        {
            if(seq_dat[curSeq].site[curSite] == '1') {
                baseCount1[(curSite+2) % 3 + 1]++;
            }
            else if(seq_dat[curSeq].site[curSite] == '2') {
                baseCount2[(curSite+2) % 3 + 1]++;
            }
            else if(seq_dat[curSeq].site[curSite] == '3') {
                baseCount3[(curSite+2) % 3 + 1]++;
            }
            else if(seq_dat[curSeq].site[curSite] == '0') {
                baseCount0[(curSite+2) % 3 + 1]++;
            }
        }
        for(i=1; i<=3; i++)
        {
            if(baseCount0[i] != seq_dat[curSeq].baseCount0[i])
            {
                printf("Error%ld in baseCount0 expected %lu found %lu for codPos %lu in seq %lu\n", pos, baseCount0[i], seq_dat[curSeq].baseCount0[i], i, curSeq);
                halt;
            }	
            if(baseCount1[i] != seq_dat[curSeq].baseCount1[i])
            {
                printf("Error%ld in baseCount1 expected %lu found %lu for codPos %lu in seq %lu\n", pos, baseCount1[i], seq_dat[curSeq].baseCount1[i], i, curSeq);
                halt;
            }
            if(baseCount2[i] != seq_dat[curSeq].baseCount2[i])
            {
                printf("Error%ld in baseCount2 expected %lu found %lu for codPos %lu in seq %lu\n", pos, baseCount2[i], seq_dat[curSeq].baseCount2[i], i, curSeq);
                halt;
            }
            if(baseCount1[i] != seq_dat[curSeq].baseCount3[i])
            {
                printf("Error%ld in baseCount3 expected %lu found %lu for codPos %lu in seq %lu\n", pos, baseCount3[i], seq_dat[curSeq].baseCount3[i], i, curSeq);
                halt;
            }
        }
    }
    return 0;
}
/*--------------------------------------------------------------------------------*/
