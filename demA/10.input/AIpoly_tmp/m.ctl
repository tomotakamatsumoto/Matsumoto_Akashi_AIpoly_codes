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
repNum						= 1		                /*	- number of replicates*/
initRunGen					= 0		                /*	- number of generations for the initial run*/
preRunGen					= 0	            	/*	- number of generations for the prerun between replicates*/
repRunGen					= 7230		            	/*	- number of generations for each replicate*/
uswitch_A					= 7220		            	/*	- number of generations after which mutation rates switch to "b" type*/
uswitch_B					= 9999999		            	/*	- number of generations after which mutation rates switch to "b" type*/
/*	- population size in burn in */
prerunSeqNum					= 1000
/*	- population size before parameter switch */
initSeqNum					= 1000		                /*	- number of sequences at the start of the simulation. if there is no population fluctuation the population remains at this size*/
/*	- population size after parameter switch */
SeqNum_a					= 100		                /*	- number of sequences at the start of the simulation. if there is no population fluctuation the population remains at this size*/
SeqNum_b					= 100		                /*	- number of sequences at the start of the simulation. if there is no population fluctuation the population remains at this size*/
seqLen						= 100000		            	/*	- length of sequences. has to be a multiple of 3*/
sampleNum					= 100		                /*	- number of sequences to be sampled from the population every replicate*/
/*	- recombination rate before parameter switch */
c							= 0.001                	/*	- the recombination rate per site per generation*/
/*	- recombination rate after parameter switch */
c_b							= 0.001                	/*	- the recombination rate per site per generation*/

/*	- basal mutation rate before parameter switch */
u_a							= 0.0000025            	/*		- the basal mutation rate per site per generation*/
/*	- basal mutation rate after parameter switch */
u_b							= 0.0000025	            	/*		- the basal mutation rate per site per generation*/

/*	- mutation rate between each nucleotide is "basal mutation rate" * "relatie mutation rate" */

/*	- relative mutation rate before parameter switch */
u01_a							= 1.000	            	/*		- the 0 to 1 mutation rate per site per generation*/
u02_a							= 1.000	            	/*		- the 0 to 2 mutation rate per site per generation*/
u03_a							= 1.000	            	/*		- the 0 to 3 mutation rate per site per generation*/
u10_a							= 1.000	            	/*		- the 1 to 0 mutation rate per site per generation*/
u12_a							= 1.000	            	/*		- the 1 to 2 mutation rate per site per generation*/
u13_a							= 1.000	            	/*		- the 1 to 3 mutation rate per site per generation*/
u20_a							= 1.000	            	/*		- the 2 to 0 mutation rate per site per generation*/
u21_a							= 1.000	            	/*		- the 2 to 1 mutation rate per site per generation*/
u23_a							= 1.000	            	/*		- the 2 to 3 mutation rate per site per generation*/
u30_a							= 1.000	            	/*		- the 3 to 0 mutation rate per site per generation*/
u31_a							= 1.000              	/*		- the 3 to 1 mutation rate per site per generation*/
u32_a							= 1.000	            	/*		- the 3 to 2 mutation rate per site per generation*/
/*	- relative mutation rate after parameter switch */
u01_b							= 1.000	            	/*		- the 0 to 1 mutation rate per site per generation*/
u02_b							= 1.000	            	/*		- the 0 to 2 mutation rate per site per generation*/
u03_b							= 1.000	            	/*		- the 0 to 3 mutation rate per site per generation*/
u10_b							= 1.000	            	/*		- the 1 to 0 mutation rate per site per generation*/
u12_b							= 1.000	            	/*		- the 1 to 2 mutation rate per site per generation*/
u13_b							= 1.000	            	/*		- the 1 to 3 mutation rate per site per generation*/
u20_b							= 1.000	            	/*		- the 2 to 0 mutation rate per site per generation*/
u21_b							= 1.000	            	/*		- the 2 to 1 mutation rate per site per generation*/
u23_b							= 1.000	            	/*		- the 2 to 3 mutation rate per site per generation*/
u30_b							= 1.000	            	/*		- the 3 to 0 mutation rate per site per generation*/
u31_b							= 1.000              	/*		- the 3 to 1 mutation rate per site per generation*/
u32_b							= 1.000	            	/*		- the 3 to 2 mutation rate per site per generation*/

/*	- fitness of each mutation type */
fitness_p							= 1.0000	            	/*		- fitness of preferred codon */
fitness_u							= 0.99957	            	/*		- fitness of unprefered codon*/

/*	- fitness of each codon in burn in */


idum_init					= -15000						/*  - initialize random number generator */

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
/*---------------------------------------------------------------------------------------*/

/*	- mutation type */
/*	- CX_Y is mutation from X to Y, numbers in X and Y are based on the codon table (e.g. 1: "TTT", 2:"TTC", ..., 64: "GGG") */
/*  - "p": silent preferred, "u": silent unpreferred, "r": replacement, "s": stop codon */

C1_1_							= "0" 
C1_2_							= "1"
C1_3_							= "0"
C1_4_							= "1"
C1_5_							= "1"
C1_6_							= "2"
C1_7_							= "1"
C1_8_							= "2"
C1_9_							= "0"
C1_10_							= "1"
C1_11_							= "0"
C1_12_							= "1"
C1_13_							= "1"
C1_14_							= "2"
C1_15_							= "1"
C1_16_							= "2"
C1_17_							= "1"
C1_18_							= "2"
C1_19_							= "1"
C1_20_							= "2"
C1_21_							= "2"
C1_22_							= "3"
C1_23_							= "2"
C1_24_							= "3"
C1_25_							= "1"
C1_26_							= "2"
C1_27_							= "1"
C1_28_							= "2"
C1_29_							= "2"
C1_30_							= "3"
C1_31_							= "2"
C1_32_							= "3"
C1_33_							= "0"
C1_34_							= "1"
C1_35_							= "0"
C1_36_							= "1"
C1_37_							= "1"
C1_38_							= "2"
C1_39_							= "1"
C1_40_							= "2"
C1_41_							= "0"
C1_42_							= "1"
C1_43_							= "0"
C1_44_							= "1"
C1_45_							= "1"
C1_46_							= "2"
C1_47_							= "1"
C1_48_							= "2"
C1_49_							= "1"
C1_50_							= "2"
C1_51_							= "1"
C1_52_							= "2"
C1_53_							= "2"
C1_54_							= "3"
C1_55_							= "2"
C1_56_							= "3"
C1_57_							= "1"
C1_58_							= "2"
C1_59_							= "1"
C1_60_							= "2"
C1_61_							= "2"
C1_62_							= "3"
C1_63_							= "2"
C1_64_							= "3"
	            	
