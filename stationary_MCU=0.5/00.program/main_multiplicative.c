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
#include "Userpref.h"
#include "header2.h"

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
	get_static_w_dat(userPrefp);
	
	/*writes the simulation information as header of output file
	also outputs the information to the screen*/
	writeFsimInfo(userPrefp);
	
	/*initrun for system to reach equilibrium,
	DO_RESET - reset data being tracked. Data being tracked include
	counts of fixations, polymorphisms, multihits*/
	printf("Initrun for %lu generation\n", userPrefp->initRunGen);
	doFsim(userPrefp, userPrefp->initRunGen, DO_RESET, INIT_RUN);

	/*initialize timer*/
	getTime(0);
	
	/*loop for repNum*/
	for(curRepNum=1; curRepNum<=userPrefp->repNum; curRepNum++)
	{
		printf("Rep %lu of %lu\n", curRepNum, userPrefp->repNum);
		
		/*prerun for independence between replicates reset counts of data being tracked after the run - DO_RESET*/
		printf("\tPreRun for %lu generations\n", userPrefp->preRunGen); 
		doFsim(userPrefp, userPrefp->preRunGen, DO_RESET, PRE_RUN);
		
		/*replicate run during which data is collected
		dont reset counts of data being tracked- DONT_RESET
		the data are reset at the end of the prerun, don't reset refers to after the replicate
		reset should not be done prior to writerepdata */ 
		//		printf("\tRepRun for %lu generations\n", repRunGen);
		doFsim(userPrefp, userPrefp->repRunGen, DONT_RESET, REP_RUN);
		
		/*collect data for the current replicate*/
		writeRepData(userPrefp);
		//		printf("\tReplicate %ld done in %lus\n", curRepNum, repTime);
	}
	
	/*free memory and close files*/
	endSim(userPrefp);
	printf("\nfinitho\n");
//	getTime(1);
	end_time = clock();
	getTime_print(userPrefp, 1);
	return 0;
}
/*--------------------------------------------------------------------------------*/
void get_static_w_dat(struct userPref_str *userPrefp)
/*
	store fitness values in array given the number of 1 states
	basically a lookup table for fitness
*/
{
	double s_val = 0;
	long i, debug = 1;
	
	w_dat_static = (double *) memAlloc (userPrefp->seqLen + 1, sizeof(double), "w_dat_static");
	s_val = userPrefp->Nes / (double) 4200;
	w_dat_static[userPrefp->seqLen] = 1.0;
	for (i = userPrefp->seqLen-1; i >= 0; i--)
	{
		w_dat_static[i] = w_dat_static[i + 1] * (1.0 - s_val); 
	}
	if (debug == 1)
	{
		printf("Nes_val: %0.3f\ns_val: %0.6f\n", userPrefp->Nes, s_val);
		printf("num1\tw_val\n");
		for (i=0; i<=userPrefp->seqLen; i++)
			printf("%0ld\t%0.6f\n", i, w_dat_static[i]);
	}
	
	return;
}

/*--------------------------------------------------------------------------------*/
/*the main alogorithm involved in the fsim: mutation, recombination, reproduction
coutfix  are done in that order for the required number of generations
genNum - the number of generations to be run
flagReset - flag which says whether to reset the counts for the current call of fsim
repType - the flag which differentiates between the three types of calls - initrun, prerun and reprun*/
long doFsim(struct userPref_str *userPrefp, long genNum, long flagReset, long repType)
{	
	
	getPopScenario(userPrefp, repType, genNum);
	
	for (curGenNum=1; curGenNum <= genNum; curGenNum++)
	{
		#ifdef DEBUGMODE
				printf("\tGeneration %lu\n", curGenNum);
				halt;
		#endif
		
//		curPopSize = popSizeArray[curGenNum];
		doMutations(userPrefp);
		
		#ifdef DEBUGMODE
				checkBaseCounts(1);
		#endif
		
		doRecombination(userPrefp);
		
		#ifdef DEBUGMODE
				checkBaseCounts(2);
		#endif
		
		getNextGenSeqs_MCU(userPrefp);			/* HAmod */
		
		updateCounts(userPrefp);
		
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
			for(curPopArrayCell=1;curPopArrayCell<=genNum;curPopArrayCell++)
			{
				popSizeArray[curPopArrayCell] = userPrefp->initSeqNum;
			}
			break;
		}
		/*PreRun*/
		case 1:
		{
			for(curPopArrayCell=1;curPopArrayCell<=genNum;curPopArrayCell++)
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
			for(curPopArrayCell=1;curPopArrayCell<=genNum;curPopArrayCell++)
			{
				popSizeArray[curPopArrayCell] = userPrefp->initSeqNum;
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
		siteData[curSite].fixNum01=0;
		siteData[curSite].fixNum10=0;
		siteData[curSite].mutNum01=0;
		siteData[curSite].mutNum10=0;
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
long doMutations(struct userPref_str *userPrefp)
{
	long curCodPos=0;							/*current codon position, 1, 2=rep, 3=sil*/	
	long curSite=0;					/*variable to loop through all flags, also used to implement mutation*/
	long curSeq=0;					/*current sequence being mutated*/
	long curHitSite=0;				/*Index of array storing hitsites*/
	double expHits10=0;						/*expected number of hits at 1 sites for a given codpos*/
	double expHits01=0;						/*expected number of hits at 0 sites for a given codpos*/
	long numHits10=0;				/*poisson deviate of the expected number for 1 sites*/
	long numHits01=0;				/*poisson deviate of the expected number for 0 sites*/
	long numChanges10=0;				/*effective number of changes from 1 to 0*/
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
	for (curSeq=1; curSeq <= userPrefp->initSeqNum; curSeq++)
	{
		/* Loop over positions 1, 2, 3 of codons */
		for (curCodPos=1; curCodPos<=3; curCodPos++)
		{	 
			expHits10 = expHits01 = 0.0;
			
			/* get expected number of hits for current cod pos */
			expHits10 = userPrefp->u10 * seq_dat[curSeq].baseCount1[curCodPos];
			expHits01 = userPrefp->u01 * seq_dat[curSeq].baseCount0[curCodPos];
			numHits10 = numHits01 = 0;
			
			/*take poisson deviates of expected numbers */
 			numHits10 = (long) poidev(expHits10);
 			numHits01 = (long) poidev(expHits01);
						
			/* limit the poisson deviates to the logical range (0 to basect for the ancestral seq)  this doesn't seem to be a problem */
			// 			if(numHits10 < 0)
			// 			{
			// 				printf("$$$  numHits10 < 0\n");
			// 				numHits10 = 0;
			// 			}
			// 			if(numHits01 < 0)
			// 			{
			// 				printf("$$$  numHits01 < 0\n");
			// 				numHits01 = 0;
			// 			}

			#ifdef DEBUGMODE_1
				printf("\n** start loop **\n");
				printf("curSeq: %0ld\tcurCodPos: %0ld\n", curSeq, curCodPos);
				printf("numHits10: %0ld\nnumHits01: %0ld\n\n", numHits10, numHits01);
			#endif

			if ((numHits10 > 0) || (numHits01 > 0))
			{
				if(numHits10 > seq_dat[curSeq].baseCount1[curCodPos])
					numHits10 = seq_dat[curSeq].baseCount1[curCodPos];
					
				/*pick sites at which to change the states*/
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
					
					#ifdef VERBOSE
						printf("10 curSeq=%lu curCodPos=%lu\n", curSeq, curCodPos);
						printf("baseCount1 = %lu baseCount0 = %lu\n", seq_dat[curSeq].baseCount1[curCodPos], seq_dat[curSeq].baseCount0[curCodPos]);
					#endif
				}
				if (numHits01 > 0)
				{
					if(numHits01 > seq_dat[curSeq].baseCount0[curCodPos])
						numHits01 = seq_dat[curSeq].baseCount0[curCodPos];
						
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
					
					#ifdef VERBOSE
						printf("01 curSeq=%lu curCodPos=%lu\n", curSeq, curCodPos);
						printf("baseCount1 = %lu baseCount0 = %lu\n", seq_dat[curSeq].baseCount1[curCodPos], seq_dat[curSeq].baseCount0[curCodPos]);
					#endif
				}
				
				/*change sites and update counts for 1->0 mutations also update siteData mutNum10*/
				if (numHits10 > 0)
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
							if (seq_dat[curSeq].site[curSite] != '1')		/* HA speed cut */
							{
								printf("mutation error 10 hit at -%c- site\n", seq_dat[curSeq].site[curSite] );
								halt;
							}
						#endif
						
						seq_dat[curSeq].site[curSite] = '0';
						siteData[curSite].mutNum10++;
						siteData[curSite].freqs1--;
						siteData[curSite].freqs0++;
						
						/* sample new sVal for the sArray for each newly arisen mutation */
						//	sArray[curSite]=sampleSval(curSite, 1);
						/* use fixed sVals fromthe sArray otherwise HAmod */
					}
					
				/*change sites and update counts for 0->1 mutations also update siteData mutNum01*/
				if (numHits01 > 0)
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
						/*get new sVal for the sArray*/
						//  sArray[curSite]=sampleSval(curSite, 0);
						/* HAmod */
					}
					
				//	printf("!!!  processing basect changes !!!\n");
				numChanges10 = 0;		
				numChanges10 = numHits10 - numHits01;
				seq_dat[curSeq].baseCount1[curCodPos] -= numChanges10;
				seq_dat[curSeq].baseCount0[curCodPos] += numChanges10;		
			}

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
				printf("Very High mutation rate / site\n");
				printf("origSite    : %c\n", origSite);
				printf("codPos      : %0ld\n", codPos);
				printf("basecount0  : %0ld\n", seq_dat[curSeq].baseCount0[codPos]);
				printf("tot_1_ct    : %0ld\n", seq_dat[curSeq].tot_1_ct);
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
		curSite = getPosNthChar_v2(userPrefp, seq_dat[curSeq].site, codPos, hitCod, origSite);
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
long getPosNthChar(struct userPref_str *userPrefp, char *site, long codPos, long hitCod, char origSite)
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
long getPosNthChar_v2(struct userPref_str *userPrefp, char *site, long codPos, long hitCod, char origSite)
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
		if (site[curSite] == origSite)
		{
			curCod++;
			if (curCod == hitCod)
				return curSite;
		}
		/*increment current site by 3 to look at codPos in the next codon*/
		curSite += 3;
	}
	printf("\ngetPosNthChar_v2 error\n");						
	printf("\tfound a total of %lu %c's in the sequence: looking for %0ldth occurrence\n", curCod, origSite, hitCod);						
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
long doRecombination(struct userPref_str *userPrefp)
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
	expCrossOvers = userPrefp->c * userPrefp->seqLen;
	
	/*The maximum number of crossovers is one less than the sequence length
	so if there are only 2 bases in the sequence the maximum number of 
	crossovers =1*/
	maxCrossOvers = userPrefp->seqLen - 1;

	/*seqIndex array is used to pair up sequences. It is initialiezed with
	numbers 1 to userPrefp->initSeqNum*/	
	for(curSeq=1; curSeq<= userPrefp->initSeqNum; curSeq++)
		seqIndex[curSeq] = curSeq;
		
	/*The indices of the pairs of sequences that are going to be recombined are
	stored at positions 1 and 2 in the seqIndex array*/	
	curSeq = 1;
	curSeqPlusOne = 2;
	
	/*repeat loop until the second of the pair is the last cell in the array*/
	while (curSeqPlusOne <= userPrefp->initSeqNum)
	{
		/*find an index to be swapped with the index at curSeqPlusOne position in the 
		seqIndex array. Randomly pick a cell from all the remaining indices*/
		nextSeq = getRandLong_(curSeqPlusOne, userPrefp->initSeqNum);
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
						if(seq_dat[curRecSeq1].site[curSite] == '1')
						{
							if (seq_dat[curRecSeq2].site[curSite] == '0')
							{
								/*site in first seq is 1 and in second is 0
								so swap and update counts*/
								curCodPos=sitePosToCodPos(curSite);
								seq_dat[curRecSeq1].site[curSite] ='0';
								seq_dat[curRecSeq1].baseCount0[curCodPos]++;
								seq_dat[curRecSeq1].baseCount1[curCodPos]--;
								seq_dat[curRecSeq2].site[curSite] = '1';
								seq_dat[curRecSeq2].baseCount0[curCodPos]--;
								seq_dat[curRecSeq2].baseCount1[curCodPos]++;
							}
						}
						else 
						{
							if (seq_dat[curRecSeq2].site[curSite] == '1')
							{
								/*site in first seq is 0 and in second is 1
								so swap and update counts*/
								curCodPos=sitePosToCodPos(curSite);
								seq_dat[curRecSeq1].site[curSite] = '1';
								seq_dat[curRecSeq1].baseCount0[curCodPos]--;
								seq_dat[curRecSeq1].baseCount1[curCodPos]++;
								seq_dat[curRecSeq2].site[curSite] = '0';
								seq_dat[curRecSeq2].baseCount0[curCodPos]++;
								seq_dat[curRecSeq2].baseCount1[curCodPos]--;
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
long getNextGenSeqs(struct userPref_str *userPrefp)
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
	
	/*Initialize weights to 1.0 each*/
	for (curSeq=1; curSeq <= userPrefp->initSeqNum; curSeq++)
	{
		wData[curSeq] = 1.0;
	}
	/*Loops through all variable sites and finds the weights for each individual*/
	for (curSite=1; curSite <= userPrefp->seqLen; curSite++)
	{	 
		if (polySite[curSite])
		{
			for (curSeq=1; curSeq<=userPrefp->initSeqNum; curSeq++)
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
	for(curSeq=1;curSeq<=userPrefp->initSeqNum;curSeq++)
	{	
		wSum += wData[curSeq];
	}
	/*expected frequencies in the nextgeneration = relative weights*/
	for(curSeq=1;curSeq<=userPrefp->initSeqNum;curSeq++)
	{	
		freqs[curSeq]=wData[curSeq]/wSum;
	}
	/*the population size for a given generation is given in the popSizeArray*/
	newPopSize = popSizeArray[curGenNum+1];
	/*take multinomial deviates on the expected frequencies*/
	multdev(freqs, userPrefp->initSeqNum, newPopSize, newFreqs);
	/*create new generation of individuals from the new frequences from 
	 multinomial sampling*/
	curSeq=1;
	for(curParent=1;curParent<=userPrefp->initSeqNum;curParent++)
	{
		for(curChild=1;curChild<=newFreqs[curParent];curChild++)
		{
//			new_seq_dat[curSeq] = seq_dat[curParent];			/* HAcheck */
			curSeq++;
		}
	}
	/*check if sampling is working correctly*/
	if ((curSeq-1) != newPopSize)
		errorOut(("Invalid number of new sequences(%lu!=%lu) in getNextGenSeqs", curSeq, newPopSize));
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
//	double wData[MAXSEQNUM+1];				/*array that stores the weight for each sequence*/
//	double freqs[MAXSEQNUM+1];				/*expected frequencies of current individuals*/
//	long newFreqs[MAXSEQNUM+1];			/*obtained frequencies of offsprings for current generation*/
long getNextGenSeqs_MCU(struct userPref_str *userPrefp)
{
	long curSeq, curCodPos;						/*refers to current sequence*/
//	long curSite, curCodPos;					
	long curParent;				/*refers to the current parent sequence*/
	long curChild;					/*refers to all child sequences of a given parent*/
//	long curSite;					/*refers to current site*/
//	long newPopSize;				/*the new population size*/
	double wSum;							/*sum of wData*/
	void *tempPV;							/*temporary pointer variable for swapping sequence pointers*/
	
	#ifdef DEBUGMODE
		printf("Finding next gen seqs\n");
		halt;
	#endif	
		
	/* get fitness for each seq and calc ave */
	wSum=0.0;
	for (curSeq=1; curSeq <= userPrefp->initSeqNum; curSeq++)
	{
//		wData[curSeq] = 0.0;
		seq_dat[curSeq].tot_1_ct = 0;
		for (curCodPos=1; curCodPos<=3; curCodPos++)
			seq_dat[curSeq].tot_1_ct += seq_dat[curSeq].baseCount1[curCodPos];
		wData[curSeq] = w_dat_static[seq_dat[curSeq].tot_1_ct];	
		wSum += wData[curSeq];			
	}
	
	/*expected frequencies in the nextgeneration = relative weights*/
	for (curSeq=1; curSeq<=userPrefp->initSeqNum; curSeq++)
		freqs[curSeq] = wData[curSeq] / wSum;

	/*the population size for a given generation is given in the popSizeArray*/
//	newPopSize = popSizeArray[curGenNum+1];
	
	/*take multinomial deviates on the expected frequencies*/
	multdev(freqs, userPrefp->initSeqNum, userPrefp->initSeqNum, newFreqs);
	
	curSeq=0;
	for (curParent=1; curParent <= userPrefp->initSeqNum; curParent++)
		for (curChild=1; curChild <= newFreqs[curParent]; curChild++)
			curSeq++;
			
	/*check if sampling is working correctly*/
	if (curSeq != userPrefp->initSeqNum)
		errorOut(("Invalid number of new sequences(%lu!=%lu) in getNextGenSeqs_MCU", curSeq, userPrefp->initSeqNum));
			
	/* option 1 - anoop's original pointer method */
	/*create new generation of individuals from the new frequences from multinomial sampling*/
	curSeq=0;
	for (curParent=1; curParent<=userPrefp->initSeqNum; curParent++)
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
long updateCounts(struct userPref_str *userPrefp)
{
	long curSite;					/*refers to the current site*/
	long curSeq;					/*refers to the current sequence*/
	long baseCount0;				/*number of 0 bases at curSite across all seqs*/
	long baseCount1;				/*number of 1 bases at curSite across all seqs*/
	
#ifdef DEBUGMODE
	printf("Updating Counts\n");
	halt;
#endif
	/*Loop across all the sites*/
	for(curSite=1; curSite <= userPrefp->seqLen; curSite++)
	{
		/*check only if the polySite flag is set*/
		if(polySite[curSite] == 1)
		{
			/*Find the frequencies of 0s and 1s at the site*/
			baseCount0=0;
			baseCount1=0;
			for (curSeq=1; curSeq<=userPrefp->initSeqNum; curSeq++)
			{
				if (seq_dat[curSeq].site[curSite] == '0')
					baseCount0++;
				else
					baseCount1++;
			}
			/*update the counts in the siteData structure*/
			siteData[curSite].freqs0 = baseCount0;
			siteData[curSite].freqs1 = baseCount1;
			/*Check if any of the basecounts is 0*/
			if(baseCount1==0)
			{
				/*if so check if the ancestral site is same as the base whose count is zero
				if so a fixation. here check if 0 has been fixed*/
				if(siteData[curSite].ancSite == '1')
				{
					/*if fixation update counts in the siteData structure*/
					siteData[curSite].ancSite = '0';
					siteData[curSite].fixNum10++;
					/*store the svalues positions and fixed site to be output in the writerepdata*/
					popFixNum++;
					fixSite[popFixNum] = '0';
					fixSitePos[popFixNum] = curSite;
					fixSiteSval[popFixNum] = sArray[curSite];
					/*at every fixation resample for the sArray*/
					sArray[curSite] = 0.0;
				}
				polySite[curSite] = 0;
			}
			else if (baseCount0 == 0)
			{
				/*if so check if the ancestral site is same as the base whose count is zero
				if so a fixation. here check if 1 has been fixed*/
				if (siteData[curSite].ancSite == '0')
				{
					/*if fixation update counts in the siteData structure*/
					siteData[curSite].ancSite = '1';
					siteData[curSite].fixNum01++;
					/*store the svalues positions and fixed site to be output in the writerepdata*/
					popFixNum++;
					fixSite[popFixNum] = '1';
					fixSitePos[popFixNum] = curSite;
					fixSiteSval[popFixNum] = sArray[curSite];
					/*at every fixation resample for the sArray*/
					sArray[curSite] = 0.0;
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
long writeRepData(struct userPref_str *userPrefp)
{
	long i, j, sum, sum01, sum10;						/*for loop*/
	long numPopPeriods;			/*number of periods for the population fluctuation scenario*/
	long startPeriod;				/*starting generation of a population period*/
	long endPeriod;				/*ending generation of a population period*/
	long curSeq;					/*refers to the current sampled sequence*/
	double ave_MCU;
	
	/*sample sequences and collect data from the samples*/
	getRepData(userPrefp);
	
//	printf("\tWriting Replicate information to file\n");
	/*Outputs the counts obtained during the current replicate*/
	fprintf(fpOut, "Replicate:%6ld\n", curRepNum);
	fprintf(fpOut, "--- seq_dat[curSeq].tot_1_ct -------------\n");
	sum = 0;
	for(curSeq=1; curSeq<= userPrefp->initSeqNum; curSeq++)
	{
//		fprintf(fpOut, "%0ld\n", seq_dat[curSeq].tot_1_ct);
		sum += seq_dat[curSeq].tot_1_ct;
	}
	ave_MCU = (double) sum / (double) userPrefp->initSeqNum;
	fprintf(fpOut, "ave: %0.4f\n", ave_MCU);
	fprintf(fpOut, "----------------\n");
	sum01 = sum10 = 0;
	sum01 += sampleFix01[1] + sampleFix01[2] + sampleFix01[3];
	sum10 += sampleFix10[1] + sampleFix10[2] + sampleFix10[3];
	fprintf(fpOut, "sampleFix01        = %lu\n", 		sum01);
	fprintf(fpOut, "sampleFix10        = %lu\n", 		sum10);
	fprintf(fpOut, "sampleFix01[1]     = %lu\n", 		sampleFix01[1]	);
	fprintf(fpOut, "sampleFix01[2]     = %lu\n", 		sampleFix01[2]	);
	fprintf(fpOut, "sampleFix01[3]     = %lu\n", 		sampleFix01[3]	);
	fprintf(fpOut, "sampleFix10[1]     = %lu\n", 		sampleFix10[1]	);
	fprintf(fpOut, "sampleFix10[2]     = %lu\n", 		sampleFix10[2]	);
	fprintf(fpOut, "sampleFix10[3]     = %lu\n", 		sampleFix10[3]	);
	sum01 = sum10 = 0;
	sum01 += sampleSeg01[1] + sampleSeg01[2] + sampleSeg01[3];
	sum10 += sampleSeg10[1] + sampleSeg10[2] + sampleSeg10[3];
	fprintf(fpOut, "sampleSeg01        = %lu\n", 		sum01);
	fprintf(fpOut, "sampleSeg10        = %lu\n", 		sum10);
	fprintf(fpOut, "sampleSeg01[1]     = %lu\n", 		sampleSeg01[1]	);
	fprintf(fpOut, "sampleSeg01[2]     = %lu\n", 		sampleSeg01[2]	);
	fprintf(fpOut, "sampleSeg01[3]     = %lu\n", 		sampleSeg01[3]	);
	fprintf(fpOut, "sampleSeg10[1]     = %lu\n", 		sampleSeg10[1]	);
	fprintf(fpOut, "sampleSeg10[2]     = %lu\n", 		sampleSeg10[2]	);
	fprintf(fpOut, "sampleSeg10[3]     = %lu\n", 		sampleSeg10[3]	);
	fprintf(fpOut, "popRecNum          = %lu\n", 		popRecNum		);
	fprintf(fpOut, "genMultiHitCount   = %lu\n", 		genMultiHitCount);
	fprintf(fpOut, "segMultiHitCount   = %lu\n", 		segMultiHitCount);
	fprintf(fpOut, "repTime            = %lu\n", 		repTime			);
	fprintf(fpOut, "\n");
	
	sum = 0;
	sum += sampleSeg01[1] + sampleSeg01[2] + sampleSeg01[3];
	sum += sampleSeg10[1] + sampleSeg10[2] + sampleSeg10[3];
	
	fprintf(fpOut2, "%0ld\t", curRepNum);
	fprintf(fpOut2, "%lu\t", sum);
	fprintf(fpOut2, "%lu\t", sampleSeg01[1]);
	fprintf(fpOut2, "%lu\t", sampleSeg01[2]);
	fprintf(fpOut2, "%lu\t", sampleSeg01[3]);
	fprintf(fpOut2, "%lu\t", sampleSeg10[1]);
	fprintf(fpOut2, "%lu\t", sampleSeg10[2]);
	fprintf(fpOut2, "%lu\t", sampleSeg10[3]);
	fprintf(fpOut2, "\n");
	
	fprintf(fpOut_segsites, "%lu\n", sum);
	fprintf(fpOut_MCU, "%0.4f\n", ave_MCU);
	
	/*outputs the population fluctuation scenario for the prerun*/
	fprintf(fpOut, "Population fluctuation for preRun(startGen:endGen:popSize)\n");
	fprintf(fpOut, "----------------------------------------------------------\n");
	/*gets the number of periods*/
	numPopPeriods=1;
	i=2;
	while(i<= userPrefp->preRunGen)
	{
		if(preRunPopFluct[i]!=preRunPopFluct[i-1])
			numPopPeriods++;
		i++;
	}
	fprintf(fpOut, "%lu\n", numPopPeriods);
	/*prints out the population fluctuation for the last period*/
	i=2;
	startPeriod=1;
	numPopPeriods=1;
	while (i <= userPrefp->preRunGen)
	{
		if(preRunPopFluct[i]!=preRunPopFluct[i-1])
		{
			endPeriod=i-1;
			fprintf(fpOut, "%lu:%lu:%lu ", startPeriod, endPeriod, preRunPopFluct[startPeriod]);
			startPeriod=i;
			numPopPeriods++;
		}
		i++;
	}
	if(numPopPeriods==1)
		fprintf(fpOut, "%d:%lu:%lu ", 1, userPrefp->preRunGen, preRunPopFluct[1]);
	fprintf(fpOut, "\n");
	fprintf(fpOut, "\n");
	
	/*outputs the population fluctuation scenario*/
	fprintf(fpOut, "Population fluctuation for repRun(startGen:endGen:popSize)\n");
	fprintf(fpOut, "----------------------------------------------------------\n");
	
	/*gets the number of periods*/
	numPopPeriods=1;
	i=2;
	while(i<=userPrefp->repRunGen)
	{
		if(popSizeArray[i]!=popSizeArray[i-1])
			numPopPeriods++;
		i++;
	}
	fprintf(fpOut, "%lu\n", numPopPeriods);
	
	/*prints out the population fluctuation for the last period*/
	i=2;
	startPeriod=1;
	numPopPeriods=1;
	while(i<= userPrefp->repRunGen)
	{
		if(popSizeArray[i]!=popSizeArray[i-1])
		{
			endPeriod=i-1;
			fprintf(fpOut, "%lu:%lu:%lu ", startPeriod, endPeriod, popSizeArray[startPeriod]);
			startPeriod=i;
			numPopPeriods++;
		}
		i++;
	}
	if(numPopPeriods==1)
		fprintf(fpOut, "%d:%ld:%ld ", 1, userPrefp->repRunGen, popSizeArray[1]);
	fprintf(fpOut, "\n");
	fprintf(fpOut, "\n");
	fprintf(fpOut, "Fixed Svals(fixedSite:fixedSval:pos)\n");
	fprintf(fpOut, "-----------------------------------\n");
	fprintf(fpOut, "%lu\n", popFixNum);
	for(i=1;i<=popFixNum;i++)
	{
		fprintf(fpOut, "%c:%.10f:%lu ", fixSite[i], fixSiteSval[i], fixSitePos[i]);
	}
	fprintf(fpOut, "\n");
	fprintf(fpOut, "\n");
	
	/*prints out the sitedata information plus sArray*/
	fprintf(fpOut, "SiteData(ancSite:sVal:freqs0:freqs1:fix01:fix10:mut01:mut10:rec)\n");
	fprintf(fpOut, "----------------------------------------------------------------\n");
	for(i=1; i<= userPrefp->seqLen; i++)
	{
		fprintf(fpOut, "%c:", siteData[i].ancSite);
		fprintf(fpOut, "%.10f:", sArray[i]);
		fprintf(fpOut, "%lu:", siteData[i].freqs0);
		fprintf(fpOut, "%lu:", siteData[i].freqs1);
		fprintf(fpOut, "%lu:", siteData[i].fixNum01);
		fprintf(fpOut, "%lu:", siteData[i].fixNum10);
		fprintf(fpOut, "%lu:", siteData[i].mutNum01);
		fprintf(fpOut, "%lu:", siteData[i].mutNum10);
		fprintf(fpOut, "%lu ", siteData[i].recNum);
	}
	fprintf(fpOut, "\n");
	fprintf(fpOut, "\n");
	fprintf(fpOut, "sampleFixPos\n");
	fprintf(fpOut, "------------\n");
	fprintf(fpOut, "%lu\n", sampleFixNum);
	for(j=1;j<=sampleFixNum;j++)
	{
		fprintf(fpOut, "%lu ", sampleFixPos[j]);
	}
	fprintf(fpOut, "\n");
	fprintf(fpOut, "\n");
	fprintf(fpOut, "sampleSegPos\n");
	fprintf(fpOut, "------------\n");
	fprintf(fpOut, "%lu\n", sampleSegNum);
	for(j=1;j<=sampleSegNum;j++)
	{
		fprintf(fpOut, "%lu ", sampleSegPos[j]);
	}
	fprintf(fpOut, "\n");
	fprintf(fpOut, "\n");
	
	/* prints out the sample sequences */
//	for(i=1;i<=sampleNum;i++)
//	{
//		curSeq=sampleSeqs[i];
//		for(j=1;j<=sampleSegNum;j++)
//		{
//			fprintf(fpOut, "%d", seq_dat[curSeq].site[sampleSegPos[j]]);
//		}
//		fprintf(fpOut, "\n");
//	}
	fprintf(fpOut, "//\n");
	fflush(fpOut);
	fflush(fpOut2);
	fflush(fpOut_segsites);
	fflush(fpOut_MCU);
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
	long curSample;				/*refers to the current sample num*/

//	printf("\tGetting Replicate Information\n");
	/*reset counts*/
	sampleFix01[1]=sampleFix01[2]=sampleFix01[3]=0;
	sampleFix10[1]=sampleFix10[2]=sampleFix10[3]=0;
	sampleFixNum=0;
	sampleSeg01[1]=sampleSeg01[2]=sampleSeg01[3]=0;
	sampleSeg10[1]=sampleSeg10[2]=sampleSeg10[3]=0;
	sampleSegNum=0;
	
	/*picks sampleNum sequences*/
	getSample(userPrefp);
	for (curSite=1; curSite<= userPrefp->seqLen; curSite++)
	{
		/*if the site is not polymorphic in the population then no need to process it*/
		if (polySite[curSite] == 1)
		{
			freq0=freq1=0;
			for(curSample=1; curSample <= userPrefp->sampleNum; curSample++)
			{
				if(seq_dat[sampleSeqs[curSample]].site[curSite] == '0')
					freq0++;
				else
					freq1++;
			}
			/*if site is polymorphic in sample*/
			if ((freq1!=0)&&(freq0!=0))
			{
				/*depending on the ancestral state update the counters*/
				if (siteData[curSite].ancSite == '1')
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
				if ((freq1==0)&&(siteData[curSite].ancSite == '1'))
				{
					sampleFix10[sitePosToCodPos(curSite)]++;
					/*keep track of fixations*/
					sampleFixNum++;
					sampleFixPos[sampleFixNum]=curSite;
				}
				/*fixation of 1 due to sampling*/
				if ((freq0==0)&&(siteData[curSite].ancSite == '0'))
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
long getSample(struct userPref_str *userPrefp)
{
	long curSample;					/*refers to the current sample*/
	long randSeq;						/*randomly selected sequence*/
	long curSeq;						/*loop through all samples*/
	
	/*repeat loop until all random sequences are found*/
	curSample=0;
	while(curSample < userPrefp->sampleNum)
	{
		/*pick a random sequence*/
		randSeq=getRandLong_(1, userPrefp->initSeqNum);
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
long checkBaseCounts(struct userPref_str *userPrefp, long pos)
{
	long i;
	long curSeq;
	long curSite;
	long baseCount0[3+1];
	long baseCount1[3+1];
	
	for (curSeq=1; curSeq<=userPrefp->initSeqNum; curSeq++)
	{
		baseCount0[1]=baseCount0[2]=baseCount0[3]=0;
		baseCount1[1]=baseCount1[2]=baseCount1[3]=0;
		for(curSite=1; curSite <= userPrefp->seqLen; curSite++)
		{
			if(seq_dat[curSeq].site[curSite] == '1')
				baseCount1[(curSite+2) % 3 + 1]++;
			else
				baseCount0[(curSite+2) %3 + 1]++;
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
		}
	}
	return 0;
}
/*--------------------------------------------------------------------------------*/
