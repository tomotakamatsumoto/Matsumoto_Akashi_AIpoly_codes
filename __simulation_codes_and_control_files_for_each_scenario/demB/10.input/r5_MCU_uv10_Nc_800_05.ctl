/*-----------------------------------------------------------------------------------------------------------------*/
/*folder and file names*/
/*---------------------------------------------------------------------------------------*/
process_choice				= 100						/* 1=codon pref test, 2=tissue MCUsp aa props, 3=tissue MCU v expr rs, 5=plot MCUsp in each tissue, see main() for others */

comments					= "1/50th_New_Mutations_Nes=20"
version 					= "035_04"		                /*	- a string of upto 10 digits for the current version of the simulation*/
machineName					= "mini_s_a"	            /* - machine on which the simulation is running*/
outFileName					= "fsim_out"				/*	- filename of output file will be generated from outFileName '_' version and extension eg fsimII_023.out*/
outFileExtn					= "txt"		                /*	- extension of output file*/

noMultihit					= 0		                	/*	- 1=multiple hits are not permitted 0=are permitted*/
repNum						= 100		                /*	- number of replicates*/
initRunGen					= 200000		                /*	- number of generations for the initial run*/
preRunGen					= 100000		            	/*	- number of generations for the prerun between replicates*/
repRunGen					= 10000		            	/*	- number of generations for each replicate*/
uswitch					= 5000		            	/*	- number of generations after which mutation rates switch to "b" type*/
initSeqNum					= 5000		                /*	- number of sequences at the start of the simulation. if there is no population fluctuation the population remains at this size*/
SeqNum_b					= 5000		                /*	- number of sequences at the start of the simulation. if there is no population fluctuation the population remains at this size*/
seqLen						= 2400		            	/*	- length of sequences. has to be a multiple of 3*/
sampleNum					= 50		                /*	- number of sequences to be sampled from the population every replicate*/
c							= 0.000002                	/*	- the recombination rate per site per generation*/
c_b							= 0.000002                	/*	- the recombination rate per site per generation*/
u10_a							= 0.000002	            	/*		- the 1 to 0 mutation rate per site per generation*/
u01_a							= 0.000002	            	/*		- the 0 to 1 mutation rate per site per generation*/
u10_b							= 0.000004	            	/*		- the 1 to 0 mutation rate per site per generation*/
u01_b							= 0.000002	            	/*		- the 0 to 1 mutation rate per site per generation*/
Nes						= 1.096	                    /* 2Nes - selection parameter for MCU model */
Nes_b						= 1.096	                    /* 2Nes - selection parameter for MCU model */
ParaA_s						= 1.000	                    /* parameter A in CFF of silent sites */
ParaA_r						= 1.000	                    /* parameter A in CFF of replacement sites */
ParaA_s_b					= 1.000	                    /* parameter A in CFF of silent sites */
ParaA_r_b					= 1.000	                    /* parameter A in CFF of replacement sites */
idum_init					= -1184						/*  - initialize random number generator */

/*CFF, a=0.1->MCU=0.8 (simulation)*/

/*---------------------------------------------------------------------------------------*/
/*these values usually doesnt have to be changed*/
/*---------------------------------------------------------------------------------------*/
rootInputFolder		: "../10.input/"				/*path+name of the root input folder*/
rootOutputFolder	: "../20.output/"				/*path+name of the root output folder*/
dirSep 				: "/"							/*the directory separator to be used*/
/*---------------------------------------------------------------------------------------*/
/*new file preferences         Not being currently used*/
/*---------------------------------------------------------------------------------------*/
fileOpenMode		: TEXT						/*the new file open mode BINARY or TEXT*/
newLine				: LF						/*the new line encode CR, LF or CRLF where CR = \r LF = \n*/
overWriteOption		: dontOverWrite				/*overWrite = files will be overwritten if already exist
												  dontOverWrite = program will exit if files already exist 
												  askUser = program will ask user as to what should be done*/
/*---------------------------------------------------------------------------------------*/
/* tissue_choice_ex_t_num_cut	= 10			/* used in gene filter for min and max ts_expr_num:  14=18 tissues, 15=22 tissues use this tissue choice for gene filtering */
/*donot change the names of the parameters as the program will be checking them to
match parameters. order of the parameters doesnt matter. both styles of c99 commenting
are supported and comments can run into multiple lines. The fields with a colon ":" as the
separator indicates values that are less frequently changed and the fields with a "=" are
values that would probably need to be changed for every new run. FileNames, foldernames and
paths are enclosed in quotes*/