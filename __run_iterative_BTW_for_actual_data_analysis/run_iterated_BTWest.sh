# required input
 # 1: fltrst
 # 2: input sequence file of BASEML (named as collapse.MFA.mfa)
 # 3: file of polymorphic sequences within an analyzed species (named as sample_seq.txt)

# 17_make_ancestral_site_probability_BTWne_method_12_mutation_categories.pl
 # conduct BTWne analysis
 # run as perl 17_make_ancestral_site_probability_BTWne_method_12_mutation_categories.pl a b c e f
 # a~d are integers which show
  # a: total number of sequences in "colapse.MFA.mfa"
  # b: number of internal nodes in the phylogeny
  # c: order of the target node to which the ancestral state will be inferred
  # informations for c and d is in fltrst file (e.g., '1    6136  AAAAAA: AAAA (1.000)' In this case, the number of internal node is four. If you are inferring the state at 3rd node of these four, you can set d = 3)  
  # e: order of the first nucleotide sequence of the analyzing species (population) in "sample_seq.txt";
  # f: order of the last nucleotide sequence of the analyzing species (population) in "sample_seq.txt";    


# 17_make_ancestral_site_probability_BTWest_method_12_mutation_categories.pl
 # conduct BTWest analysis
 # run as perl 17_make_ancestral_site_probability_BTWne_method_12_mutation_categories.pl a b c d e f mut_TC ~ mut_GA
 # a~c, e and f are same as above and d decides the number of mutation categories
 # if you set d = 12, the SFS for the weighting will be calculated for each of 12 mutations.
 # if you set d = 3, you can specify mutations in each of three categories and then, the SFS will be calculated by summing up all mutations in each category.
 # mut_TC~mut_GA specify the mutation category of each mutation.
 
 # For example) if you set d = 3, you can assign 1, 2 or 3 for each mutation and mutations in the same category will be summed up to calculate the SFS for the weighting
 # For example) d=3;
 # For example) mut_TC=1;
 # For example) mut_TA=3;
 # For example) mut_TG=1;
 # For example) mut_CT=2;
 # For example) mut_CA=2;
 # For example) mut_CG=3;
 # For example) mut_AT=3;
 # For example) mut_AC=1;
 # For example) mut_AG=1;
 # For example) mut_GT=2;
 # For example) mut_GC=3;
 # For example) mut_GA=2;
 
 # For example) if you set d = 12, you can assign 1, 2, ..., 12 for each mutation and SFS for the weighting will be calculated for each of 12 mutations
 # For example) d=12;
 # For example) mut_TC=1;
 # For example) mut_TA=2;
 # For example) mut_TG=3;
 # For example) mut_CT=4;
 # For example) mut_CA=5;
 # For example) mut_CG=6;
 # For example) mut_AT=7;
 # For example) mut_AC=8;
 # For example) mut_AG=9;
 # For example) mut_GT=10;
 # For example) mut_GC=11;
 # For example) mut_GA=12;
 
 
# 17_make_estimated_SFS_BTW_method.pl
 # make the estimated SFS using the inferred ancestral states
 # run as perl 17_make_estimated_SFS_BTW_method.pl e f

# iteration: number of iteration in the BTW analysis (if you set iteration = 6, first analysis is done under BTWne and then BTWest analysis will be iterated 5 times).


a=6;
b=4;
c=4;
d=12;
e=15;
f=35;

mut_TC=1;
mut_TA=2;
mut_TG=3;
mut_CT=4;
mut_CA=5;
mut_CG=6;
mut_AT=7;
mut_AC=8;
mut_AG=9;
mut_GT=10;
mut_GC=11;
mut_GA=12;

iteration=2;


dirpath="EXAMPLE"     # input folder (same directory as .sh file)

cp $dirpath/fltrst fltrst;
cp $dirpath/collapse.MFA.mfa collapse.MFA.mfa;
cp $dirpath/sample_seq.txt sample_seq.txt;

perl 17_make_ancestral_site_probability_BTWne_method_12_mutation_categories.pl $a $b $c $e $f;
perl 17_make_estimated_SFS_BTW_method.pl $e $f;

rm anc_site_probs.txt;
#cp estimated_frequency_spectrum.txt estimated_frequency_spectrum_BTWne.txt;

Y=1;
while [ $Y -ne $iteration ]

do
 
 rm AWP_substitution_for_post_weighted_collapse.txt;
 
 perl 17_make_ancestral_site_probability_BTWest_method_12_mutation_categories.pl $a $b $c $d $e $f $mut_TC $mut_TA $mut_TG $mut_CT $mut_CA $mut_CG $mut_AT $mut_AC $mut_AG $mut_GT $mut_GC $mut_GA;
 rm estimated_frequency_spectrum.txt;
 
 
 perl 17_make_estimated_SFS_BTW_method.pl $e $f;
 
 rm anc_site_probs.txt;
 #cp estimated_frequency_spectrum.txt estimated_frequency_spectrum_BTWest_"$Y".txt;
 
 $((Y++));

done


rm collapse.MFA.mfa;
rm fltrst;
rm sample_seq.txt;


