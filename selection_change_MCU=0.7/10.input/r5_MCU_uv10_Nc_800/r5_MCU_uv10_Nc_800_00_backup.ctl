/*-----------------------------------------------------------------------------------------------------------------*/
/*folder and file names*/
/*---------------------------------------------------------------------------------------*/
process_choice				= 100						/* 1=codon pref test, 2=tissue MCUsp aa props, 3=tissue MCU v expr rs, 5=plot MCUsp in each tissue, see main() for others */

comments					= "1/50th_New_Mutations_Nes=20"
version 					= "035_04"		                /*	- a string of upto 10 digits for the current version of the simulation*/
machineName					= "mini_s_a"	            /* - machine on which the simulation is running*/
outFileName					= "fsim_out"				/*	- filename of output file will be generated from outFileName '_' version and extension eg fsimII_023.out*/
outFileExtn					= "txt"		                /*	- extension of output file*/

noMultihit					= 1		                	/*	- 1=multiple hits are not permitted 0=are permitted*/
repNum						= 1		                /*	- number of replicates*/
initRunGen					= 50000		                /*	- number of generations for the initial run*/
preRunGen					= 1	            	/*	- number of generations for the prerun between replicates*/
repRunGen					= 20000		            	/*	- number of generations for each replicate*/
uswitch					= 15000		            	/*	- number of generations after which mutation rates switch to "b" type*/
/*	- population size in burn in */
prerunSeqNum					= 1000
/*	- population size before parameter switch */
initSeqNum					= 1000		                /*	- number of sequences at the start of the simulation. if there is no population fluctuation the population remains at this size*/
/*	- population size after parameter switch */
SeqNum_b					= 1000		                /*	- number of sequences at the start of the simulation. if there is no population fluctuation the population remains at this size*/
seqLen						= 3000000		            	/*	- length of sequences. has to be a multiple of 3*/
sampleNum					= 50		                /*	- number of sequences to be sampled from the population every replicate*/
/*	- recombination rate before parameter switch */
c							= 0.01                	/*	- the recombination rate per site per generation*/
/*	- recombination rate after parameter switch */
c_b							= 0.01                	/*	- the recombination rate per site per generation*/

/*	- basal mutation rate before parameter switch */
u_a							= 0.00001	            	/*		- the basal mutation rate per site per generation*/
/*	- basal mutation rate after parameter switch */
u_b							= 0.00001	            	/*		- the basal mutation rate per site per generation*/


/*	- mutation rate between each nucleotide is "basal mutation rate" * "relatie mutation rate" */

/*	- relative mutation rate before parameter switch */
u01_a							= 0.330	            	/*		- the 0 to 1 mutation rate per site per generation*/
u02_a							= 0.330	            	/*		- the 0 to 2 mutation rate per site per generation*/
u03_a							= 0.330	            	/*		- the 0 to 3 mutation rate per site per generation*/
u10_a							= 0.330	            	/*		- the 1 to 0 mutation rate per site per generation*/
u12_a							= 0.330	            	/*		- the 1 to 2 mutation rate per site per generation*/
u13_a							= 0.330	            	/*		- the 1 to 3 mutation rate per site per generation*/
u20_a							= 0.330	            	/*		- the 2 to 0 mutation rate per site per generation*/
u21_a							= 0.330	            	/*		- the 2 to 1 mutation rate per site per generation*/
u23_a							= 0.330	            	/*		- the 2 to 3 mutation rate per site per generation*/
u30_a							= 0.330	            	/*		- the 3 to 0 mutation rate per site per generation*/
u31_a							= 0.330              	/*		- the 3 to 1 mutation rate per site per generation*/
u32_a							= 0.330	            	/*		- the 3 to 2 mutation rate per site per generation*/
/*	- relative mutation rate after parameter switch */
u01_b							= 0.330	            	/*		- the 0 to 1 mutation rate per site per generation*/
u02_b							= 0.330            	    /*		- the 0 to 2 mutation rate per site per generation*/
u03_b							= 0.330	            	/*		- the 0 to 3 mutation rate per site per generation*/
u10_b							= 0.330	            	/*		- the 1 to 0 mutation rate per site per generation*/
u12_b							= 0.330	            	/*		- the 1 to 2 mutation rate per site per generation*/
u13_b							= 0.330	            	/*		- the 1 to 3 mutation rate per site per generation*/
u20_b							= 0.330	            	/*		- the 2 to 0 mutation rate per site per generation*/
u21_b							= 0.330	            	/*		- the 2 to 1 mutation rate per site per generation*/
u23_b							= 0.330	            	/*		- the 2 to 3 mutation rate per site per generation*/
u30_b							= 0.330	            	/*		- the 3 to 0 mutation rate per site per generation*/
u31_b							= 0.330	            	/*		- the 3 to 1 mutation rate per site per generation*/
u32_b							= 0.330	            	/*		- the 3 to 2 mutation rate per site per generation*/
/*	- fitness of each codon in burn in */
f000_a							= 1.000000	            	/*		- fitness of one TTT*/
f001_a							= 1.000000	            	/*		- fitness of one TTC*/
f002_a							= 1.000000	            	/*		- fitness of one TTA*/
f003_a							= 1.000000	            	/*		- fitness of one TTG*/
f010_a							= 1.000000	            	/*		- fitness of one TCT*/
f011_a							= 1.000000	            	/*		- fitness of one TCC*/
f012_a							= 1.000000	            	/*		- fitness of one TCA*/
f013_a							= 1.000000	            	/*		- fitness of one TCG*/
f020_a							= 1.000000	            	/*		- fitness of one TAT*/
f021_a							= 1.000000	            	/*		- fitness of one TAC*/
f022_a							= 1.000000	            	/*		- fitness of one TAA*/
f023_a							= 1.000000	            	/*		- fitness of one TAG*/
f030_a							= 1.000000	            	/*		- fitness of one TGT*/
f031_a							= 1.000000	            	/*		- fitness of one TGC*/
f032_a							= 1.000000	            	/*		- fitness of one TGA*/
f033_a							= 1.000000	            	/*		- fitness of one TGG*/
f100_a							= 1.000000	            	/*		- fitness of one CTT*/
f101_a							= 1.000000	            	/*		- fitness of one CTC*/
f102_a							= 1.000000	            	/*		- fitness of one CTA*/
f103_a							= 1.000000	            	/*		- fitness of one CTG*/
f110_a							= 1.000000	            	/*		- fitness of one CCT*/
f111_a							= 1.000000	            	/*		- fitness of one CCC*/
f112_a							= 1.000000	            	/*		- fitness of one CCA*/
f113_a							= 1.000000	            	/*		- fitness of one CCG*/
f120_a							= 1.000000	            	/*		- fitness of one CAT*/
f121_a							= 1.000000	            	/*		- fitness of one CAC*/
f122_a							= 1.000000	            	/*		- fitness of one CAA*/
f123_a							= 1.000000	            	/*		- fitness of one CAG*/
f130_a							= 1.000000	            	/*		- fitness of one CGT*/
f131_a							= 1.000000	            	/*		- fitness of one CGC*/
f132_a							= 1.000000	            	/*		- fitness of one CGA*/
f133_a							= 1.000000	            	/*		- fitness of one CGG*/
f200_a							= 1.000000	            	/*		- fitness of one ATT*/
f201_a							= 1.000000	            	/*		- fitness of one ATC*/
f202_a							= 1.000000	            	/*		- fitness of one ATA*/
f203_a							= 1.000000	            	/*		- fitness of one ATG*/
f210_a							= 1.000000	            	/*		- fitness of one ACT*/
f211_a							= 1.000000	            	/*		- fitness of one ACC*/
f212_a							= 1.000000	            	/*		- fitness of one ACA*/
f213_a							= 1.000000	            	/*		- fitness of one ACG*/
f220_a							= 1.000000	            	/*		- fitness of one AAT*/
f221_a							= 1.000000	            	/*		- fitness of one AAC*/
f222_a							= 1.000000	            	/*		- fitness of one AAA*/
f223_a							= 1.000000	            	/*		- fitness of one AAG*/
f230_a							= 1.000000	            	/*		- fitness of one AGT*/
f231_a							= 1.000000	            	/*		- fitness of one AGC*/
f232_a							= 1.000000	            	/*		- fitness of one AGA*/
f233_a							= 1.000000	            	/*		- fitness of one AGG*/
f300_a							= 1.000000	            	/*		- fitness of one GTT*/
f301_a							= 1.000000	            	/*		- fitness of one GTC*/
f302_a							= 1.000000	            	/*		- fitness of one GTA*/
f303_a							= 1.000000	            	/*		- fitness of one GTG*/
f310_a							= 1.000000	            	/*		- fitness of one GCT*/
f311_a							= 1.000000	            	/*		- fitness of one GCC*/
f312_a							= 1.000000	            	/*		- fitness of one GCA*/
f313_a							= 1.000000	            	/*		- fitness of one GCG*/
f320_a							= 1.000000	            	/*		- fitness of one GAT*/
f321_a							= 1.000000	            	/*		- fitness of one GAC*/
f322_a							= 1.000000	            	/*		- fitness of one GAA*/
f323_a							= 1.000000	            	/*		- fitness of one GAG*/
f330_a							= 1.000000	            	/*		- fitness of one GGT*/
f331_a							= 1.000000	            	/*		- fitness of one GGC*/
f332_a							= 1.000000	            	/*		- fitness of one GGA*/
f333_a							= 1.000000	            	/*		- fitness of one GGG*/

/*	- fitness of each codon before parameter switch */
f000_b							= 1.000000	            	/*		- fitness of one TTT*/
f001_b							= 0.950000	            	/*		- fitness of one TTC*/
f002_b							= 0.950000	            	/*		- fitness of one TTA*/
f003_b							= 0.950000	            	/*		- fitness of one TTG*/
f010_b							= 0.950000	            	/*		- fitness of one TCT*/
f011_b							= 0.950000	            	/*		- fitness of one TCC*/
f012_b							= 0.950000	            	/*		- fitness of one TCA*/
f013_b							= 0.950000	            	/*		- fitness of one TCG*/
f020_b							= 0.950000	            	/*		- fitness of one TAT*/
f021_b							= 0.950000	            	/*		- fitness of one TAC*/
f022_b							= 0.950000	            	/*		- fitness of one TAA*/
f023_b							= 0.950000	            	/*		- fitness of one TAG*/
f030_b							= 0.950000	            	/*		- fitness of one TGT*/
f031_b							= 0.950000	            	/*		- fitness of one TGC*/
f032_b							= 0.950000	            	/*		- fitness of one TGA*/
f033_b							= 0.950000	            	/*		- fitness of one TGG*/
f100_b							= 0.950000	            	/*		- fitness of one CTT*/
f101_b							= 0.950000	            	/*		- fitness of one CTC*/
f102_b							= 0.950000	            	/*		- fitness of one CTA*/
f103_b							= 0.950000	            	/*		- fitness of one CTG*/
f110_b							= 0.950000	            	/*		- fitness of one CCT*/
f111_b							= 1.000000	            	/*		- fitness of one CCC*/
f112_b							= 0.950000	            	/*		- fitness of one CCA*/
f113_b							= 0.950000	            	/*		- fitness of one CCG*/
f120_b							= 0.950000	            	/*		- fitness of one CAT*/
f121_b							= 0.950000	            	/*		- fitness of one CAC*/
f122_b							= 0.950000	            	/*		- fitness of one CAA*/
f123_b							= 0.950000	            	/*		- fitness of one CAG*/
f130_b							= 0.950000	            	/*		- fitness of one CGT*/
f131_b							= 0.950000	            	/*		- fitness of one CGC*/
f132_b							= 0.950000	            	/*		- fitness of one CGA*/
f133_b							= 0.950000	            	/*		- fitness of one CGG*/
f200_b							= 0.950000	            	/*		- fitness of one ATT*/
f201_b							= 0.950000	            	/*		- fitness of one ATC*/
f202_b							= 0.950000	            	/*		- fitness of one ATA*/
f203_b							= 0.950000	            	/*		- fitness of one ATG*/
f210_b							= 0.950000	            	/*		- fitness of one ACT*/
f211_b							= 0.950000	            	/*		- fitness of one ACC*/
f212_b							= 0.950000	            	/*		- fitness of one ACA*/
f213_b							= 0.950000	            	/*		- fitness of one ACG*/
f220_b							= 0.950000	            	/*		- fitness of one AAT*/
f221_b							= 0.950000	            	/*		- fitness of one AAC*/
f222_b							= 1.000000	            	/*		- fitness of one AAA*/
f223_b							= 0.950000	            	/*		- fitness of one AAG*/
f230_b							= 0.950000	            	/*		- fitness of one AGT*/
f231_b							= 0.950000	            	/*		- fitness of one AGC*/
f232_b							= 0.950000	            	/*		- fitness of one AGA*/
f233_b							= 0.950000	            	/*		- fitness of one AGG*/
f300_b							= 0.950000	            	/*		- fitness of one GTT*/
f301_b							= 0.950000	            	/*		- fitness of one GTC*/
f302_b							= 0.950000	            	/*		- fitness of one GTA*/
f303_b							= 0.950000	            	/*		- fitness of one GTG*/
f310_b							= 0.950000	            	/*		- fitness of one GCT*/
f311_b							= 0.950000	            	/*		- fitness of one GCC*/
f312_b							= 0.950000	            	/*		- fitness of one GCA*/
f313_b							= 0.950000	            	/*		- fitness of one GCG*/
f320_b							= 0.950000	            	/*		- fitness of one GAT*/
f321_b							= 0.950000	            	/*		- fitness of one GAC*/
f322_b							= 0.950000	            	/*		- fitness of one GAA*/
f323_b							= 0.950000	            	/*		- fitness of one GAG*/
f330_b							= 0.950000	            	/*		- fitness of one GGT*/
f331_b							= 0.950000	            	/*		- fitness of one GGC*/
f332_b							= 0.950000	            	/*		- fitness of one GGA*/
f333_b							= 1.000000	            	/*		- fitness of one GGG*/

/*	- fitness of each codon after parameter switch */
f000_c							= 1.000000	            	/*		- fitness of one TTT*/
f001_c							= 0.950000	            	/*		- fitness of one TTC*/
f002_c							= 0.950000	            	/*		- fitness of one TTA*/
f003_c							= 0.950000	            	/*		- fitness of one TTG*/
f010_c							= 0.950000	            	/*		- fitness of one TCT*/
f011_c							= 0.950000	            	/*		- fitness of one TCC*/
f012_c							= 0.950000	            	/*		- fitness of one TCA*/
f013_c							= 0.950000	            	/*		- fitness of one TCG*/
f020_c							= 0.950000	            	/*		- fitness of one TAT*/
f021_c							= 0.950000	            	/*		- fitness of one TAC*/
f022_c							= 0.950000	            	/*		- fitness of one TAA*/
f023_c							= 0.950000	            	/*		- fitness of one TAG*/
f030_c							= 0.950000	            	/*		- fitness of one TGT*/
f031_c							= 0.950000	            	/*		- fitness of one TGC*/
f032_c							= 0.950000	            	/*		- fitness of one TGA*/
f033_c							= 0.950000	            	/*		- fitness of one TGG*/
f100_c							= 0.950000	            	/*		- fitness of one CTT*/
f101_c							= 0.950000	            	/*		- fitness of one CTC*/
f102_c							= 0.950000	            	/*		- fitness of one CTA*/
f103_c							= 0.950000	            	/*		- fitness of one CTG*/
f110_c							= 0.950000	            	/*		- fitness of one CCT*/
f111_c							= 1.000000	            	/*		- fitness of one CCC*/
f112_c							= 0.950000	            	/*		- fitness of one CCA*/
f113_c							= 0.950000	            	/*		- fitness of one CCG*/
f120_c							= 0.950000	            	/*		- fitness of one CAT*/
f121_c							= 0.950000	            	/*		- fitness of one CAC*/
f122_c							= 0.950000	            	/*		- fitness of one CAA*/
f123_c							= 0.950000	            	/*		- fitness of one CAG*/
f130_c							= 0.950000	            	/*		- fitness of one CGT*/
f131_c							= 0.950000	            	/*		- fitness of one CGC*/
f132_c							= 0.950000	            	/*		- fitness of one CGA*/
f133_c							= 0.950000	            	/*		- fitness of one CGG*/
f200_c							= 0.950000	            	/*		- fitness of one ATT*/
f201_c							= 0.950000	            	/*		- fitness of one ATC*/
f202_c							= 0.950000	            	/*		- fitness of one ATA*/
f203_c							= 0.950000	            	/*		- fitness of one ATG*/
f210_c							= 0.950000	            	/*		- fitness of one ACT*/
f211_c							= 0.950000	            	/*		- fitness of one ACC*/
f212_c							= 0.950000	            	/*		- fitness of one ACA*/
f213_c							= 0.950000	            	/*		- fitness of one ACG*/
f220_c							= 0.950000	            	/*		- fitness of one AAT*/
f221_c							= 0.950000	            	/*		- fitness of one AAC*/
f222_c							= 1.000000	            	/*		- fitness of one AAA*/
f223_c							= 0.950000	            	/*		- fitness of one AAG*/
f230_c							= 0.950000	            	/*		- fitness of one AGT*/
f231_c							= 0.950000	            	/*		- fitness of one AGC*/
f232_c							= 0.950000	            	/*		- fitness of one AGA*/
f233_c							= 0.950000	            	/*		- fitness of one AGG*/
f300_c							= 0.950000	            	/*		- fitness of one GTT*/
f301_c							= 0.950000	            	/*		- fitness of one GTC*/
f302_c							= 0.950000	            	/*		- fitness of one GTA*/
f303_c							= 0.950000	            	/*		- fitness of one GTG*/
f310_c							= 0.950000	            	/*		- fitness of one GCT*/
f311_c							= 0.950000	            	/*		- fitness of one GCC*/
f312_c							= 0.950000	            	/*		- fitness of one GCA*/
f313_c							= 0.950000	            	/*		- fitness of one GCG*/
f320_c							= 0.950000	            	/*		- fitness of one GAT*/
f321_c							= 0.950000	            	/*		- fitness of one GAC*/
f322_c							= 0.950000	            	/*		- fitness of one GAA*/
f323_c							= 0.950000	            	/*		- fitness of one GAG*/
f330_c							= 0.950000	            	/*		- fitness of one GGT*/
f331_c							= 0.950000	            	/*		- fitness of one GGC*/
f332_c							= 0.950000	            	/*		- fitness of one GGA*/
f333_c							= 1.000000	            	/*		- fitness of one GGG*/

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