# required input
 # 1: fltrst
 # 2: input sequence file of BASEML (named as collapse.MFA.mfa)
 # 3: file of polymorphic sequences within an analyzed species (named as sample_seq.txt)

# 17_make_ancestral_site_probability_BTWne_method_12_mutation_categories.pl
 # conduct BTWne analysis
 # run as perl 17_make_ancestral_site_probability_BTWne_method_12_mutation_categories.pl a b c d
 # a~d are integers which show
  # a: total number of sequences in "colapse.MFA.mfa"
  # b: number of polymorphic sequences in "sample_seq.txt"
  # c: number of internal nodes in the phylogeny
  # d: order of the target node to which the ancestral state will be inferred
  # informations for c and d is in fltrst file (e.g., '1    6136  AAAAAA: AAAA (1.000)' In this case, the number of internal node is four. If you are inferring the state at 3rd node of these four, you can set d = 3)  
  


# 17_make_ancestral_site_probability_BTWest_method_12_mutation_categories.pl
 # conduct BTWest analysis
 # run as perl 17_make_ancestral_site_probability_BTWne_method_12_mutation_categories.pl a b c d e
 # a~d are same as above and e decides the number of mutation categories
 # if you set e = 12, the SFS for the weighting will be calculated for each of 12 mutations.
 # if you set e = 3, you can specify mutations in each of three categories and then, the SFS will be calculated by summing up all mutations in each category.
 # the code will ask the mutations in each category. so please type the number corresponds to the mutations as below
 # 1:TC
 # 2:TA
 # 3:TG
 # 4:CT
 # 5:CA
 # 6:CG
 # 7:AT
 # 8:AC
 # 9:AG
 # 10:GT
 # 11:GC
 # 12:GA
 # For example) if you set e = 3, the code will ask "enter the mutations in category 1" and then, you can type 1 3 8 9 (category 1 is AT->GC) 
 # For example)                                     "enter the mutations in category 2" and then, you can type 4 5 10 12 (category 2 is GC->AT) 
 # For example)                                     "enter the mutations in category 3" and then, you can type 2 6 7 11 (category 3 is GC->GC and AT->AT) 
 
 
# 17_make_estimated_SFS_BTW_method.pl
 # make the estimated SFS using the inferred ancestral states
 # run as perl 17_make_estimated_SFS_BTW_method.pl b
 # b is the number of polymorphic sequences as explained above




a=12;
b=10;
c=10;
d=3;
e=3;


dirpath="EXAMPLE"     # input folder (same directory as .sh file)

cp $dirpath/fltrst fltrst;
cp $dirpath/collapse.MFA.mfa collapse.MFA.mfa;
cp $dirpath/sample_seq.txt sample_seq.txt;

perl 17_make_ancestral_site_probability_BTWne_method_12_mutation_categories.pl $a $b $c $d;
perl 17_make_estimated_SFS_BTW_method.pl $b;

rm anc_site_probs.txt;
#cp estimated_frequency_spectrum.txt estimated_frequency_spectrum_BTWne.txt;

Y=1;
while [ $Y -ne 6 ]

do
 
 rm AWP_substitution_for_post_weighted_collapse.txt;
 
 perl 17_make_ancestral_site_probability_BTWest_method_12_mutation_categories.pl $a $b $c $d $e;
 rm estimated_frequency_spectrum.txt;
 
 
 perl 17_make_estimated_SFS_BTW_method.pl $b;
 
 rm anc_site_probs.txt;
 #cp estimated_frequency_spectrum.txt estimated_frequency_spectrum_BTWest_"$Y".txt;
 
 $((Y++));

done

 cat estimated_frequency_spectrum.txt  >> estimated_frequency_spectrum_iterative_BTWest_bin"$Z".txt;
 cat AWP_substitution_for_post_weighted_collapse.txt  >> AWP_substitution_iterative_BTWest_bin"$Z".txt;
 
 rm estimated_frequency_spectrum.txt;
 rm AWP_substitution_for_post_weighted_collapse.txt;

 rm collapse.MFA.mfa;
 rm fltrst;
 rm anc_site_prob_m_0%_ms.txt;
 rm site_pattern_with_higher_than_0%_ms_node.txt;
 rm site_pattern_with_lower_than_0%_ms_node.txt;

