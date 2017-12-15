1/50th_New_Mutations_Nes=20
035		                /*version	- a string of upto 10 digits for the current version of the simulation*/
04		                /*fsimIIversion - a string of upto 10 digits for the current version of fsim*/
mini_s_a	            /*machineName - machine on which the simulation is running*/
../20.output/fsim_out	/*outFileName	- filename of output file will be generated from outFileName '_' version and extension eg fsimII_023.out*/
txt		                /*outFileExtn	- extension of output file*/
1		                /*noMultihit	- 1=multiple hits are not permitted 0=are permitted*/
100		                /*repNum	- number of replicates*/
500000		            /*initRunGen	- number of generations for the initial run*/
20000		            /*preRunGen	- number of generations for the prerun between replicates*/
2000		            /*repRunGen	- number of generations for each replicate*/
500		                /*initSeqNum	- number of sequences at the start of the simulation. if there is no population fluctuation the population remains at this size*/
1200		            /*seqLen	- length of sequences. has to be a multiple of 3*/
25		                /*sampleNum	- number of sequences to be sampled from the population every replicate*/
0.001	                /*c		- the recombination rate per site per generation*/
0.0000125	            /*u10		- the 1 to 0 mutation rate per site per generation*/
0.0000125	            /*u01		- the 0 to 1 mutation rate per site per generation*/
1.0	                    /* 2Nes		- selection parameter for MCU model */
-4397					/* idum_init - initialize random number generator */

/*The first line is a comment which will be read and output 
into the fsim output file. dont add any comments to the end of the first line
whatever comes before the first end of line is assumed to be a comment
The order of the control file items is significant. DONT CHANGE ORDER. if you
want to change order then change the getCtlFile function as well*/



	readDataWithComments(fpCtl, "%s", version);
	readDataWithComments(fpCtl, "%s", fsimVersion);
	readDataWithComments(fpCtl, "%s", machineName);
	readDataWithComments(fpCtl, "%s", outFileName);
	readDataWithComments(fpCtl, "%s", outFileExtn);
	readDataWithComments(fpCtl, "%ld", &noMultiHit);
	readDataWithComments(fpCtl, "%ld", &repNum);
	readDataWithComments(fpCtl, "%ld", &initRunGen);
//	printf("initRunGen: %0ld\n", initRunGen);
	readDataWithComments(fpCtl, "%ld", &preRunGen);
//	printf("preRunGen: %0ld\n", preRunGen);
	readDataWithComments(fpCtl, "%ld", &repRunGen);
//	printf("repRunGen: %0ld\n", repRunGen);
	readDataWithComments(fpCtl, "%ld", &initSeqNum);
	readDataWithComments(fpCtl, "%ld", &seqLen);
	readDataWithComments(fpCtl, "%ld", &sampleNum);
	readDataWithComments(fpCtl, "%lf", &c);
	readDataWithComments(fpCtl, "%lf", &u10);
	readDataWithComments(fpCtl, "%lf", &u01);
	readDataWithComments(fpCtl, "%lf", &Nes2);
	readDataWithComments(fpCtl, "%ld", &idum_init);
