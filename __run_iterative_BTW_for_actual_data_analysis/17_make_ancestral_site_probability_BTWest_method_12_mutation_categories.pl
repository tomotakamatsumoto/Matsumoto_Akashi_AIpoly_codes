open (IN1,"collapse.MFA.mfa");
@str1 = <IN1>;
open (IN2,"fltrst");
@str2 = <IN2>;
open (IN3,"sample_seq.txt");
@str3 = <IN3>;

open (IN4,"estimated_frequency_spectrum.txt");
@str4 = <IN4>;

# read integers : total number of sequences (collapse seqs + outgroups), number of alleles in the focused species, number of internal nodes in the phylogeny, order of the target node to which the ancestral state is inferred
$total_seq_num = $ARGV[0];
$num_internal_nodes = $ARGV[1];
$target_node = $ARGV[2];
$sample_start = $ARGV[4];
$sample_end = $ARGV[5];

$sample_num = $sample_end - $sample_start + 1;


$mutation_category_num = $ARGV[3];

@TC = split /\t/, @str4[0];
@TA = split /\t/, @str4[1];
@TG = split /\t/, @str4[2];
@CT = split /\t/, @str4[3];
@CA = split /\t/, @str4[4];
@CG = split /\t/, @str4[5];
@AT = split /\t/, @str4[6];
@AC = split /\t/, @str4[7];
@AG = split /\t/, @str4[8];
@GT = split /\t/, @str4[9];
@GC = split /\t/, @str4[10];
@GA = split /\t/, @str4[11];

$num_TC=0;
$num_TA=0;
$num_TG=0;
$num_CT=0;
$num_CA=0;
$num_CG=0;
$num_AT=0;
$num_AC=0;
$num_AG=0;
$num_GT=0;
$num_GC=0;
$num_GA=0;

# sum the numbers of polymorphic mutations in one mutation category
for ($n1=1;$n1<=$mutation_category_num;$n1++) {
    for ($n2=1;$n2<=$sample_num;$n2++) {
        @{'num_category'.$n1}[$n2] = 0;
    }
    #$check_catenum=1;
    #$cate_1 = $ARGV[6];
    #for ($n3=7;$n3=17;$n++) {
    #    if ($ARGV[$n3] !=
    #}
    
    
    $cate_check_TC = $ARGV[6];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[6]}[$n2] = @{'num_category'.$ARGV[6]}[$n2]+@TC[$n2];
    }
    $cate_check_TA = $ARGV[7];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[7]}[$n2] = @{'num_category'.$ARGV[7]}[$n2]+@TA[$n2];
    }
    $cate_check_TG = $ARGV[8];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[8]}[$n2] = @{'num_category'.$ARGV[8]}[$n2]+@TG[$n2];
    }
    $cate_check_CT = $ARGV[9];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[9]}[$n2] = @{'num_category'.$ARGV[9]}[$n2]+@CT[$n2];
    }
    $cate_check_CA = $ARGV[10];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[10]}[$n2] = @{'num_category'.$ARGV[10]}[$n2]+@CA[$n2];
    }
    $cate_check_CG = $ARGV[11];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[11]}[$n2] = @{'num_category'.$ARGV[11]}[$n2]+@CG[$n2];
    }
    $cate_check_AT = $ARGV[12];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[12]}[$n2] = @{'num_category'.$ARGV[12]}[$n2]+@AT[$n2];
    }
    $cate_check_AC = $ARGV[13];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[13]}[$n2] = @{'num_category'.$ARGV[13]}[$n2]+@AC[$n2];
    }
    $cate_check_AG = $ARGV[14];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[14]}[$n2] = @{'num_category'.$ARGV[14]}[$n2]+@AG[$n2];
    }
    $cate_check_GT = $ARGV[15];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[15]}[$n2] = @{'num_category'.$ARGV[15]}[$n2]+@GT[$n2];
    }
    $cate_check_GC = $ARGV[16];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[16]}[$n2] = @{'num_category'.$ARGV[16]}[$n2]+@GC[$n2];
    }
    $cate_check_GA = $ARGV[17];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[17]}[$n2] = @{'num_category'.$ARGV[17]}[$n2]+@GA[$n2];
    }
}

printf "TC is in category $cate_check_TC\n";
printf "TA is in category $cate_check_TA\n";
printf "TG is in category $cate_check_TG\n";
printf "CT is in category $cate_check_CT\n";
printf "CA is in category $cate_check_CA\n";
printf "CG is in category $cate_check_CG\n";
printf "AT is in category $cate_check_AT\n";
printf "AC is in category $cate_check_AC\n";
printf "AG is in category $cate_check_AG\n";
printf "GT is in category $cate_check_GT\n";
printf "GC is in category $cate_check_GC\n";
printf "GA is in category $cate_check_GA\n";

# calculate the SFSest of each mutation category
for ($n1=1;$n1<=$mutation_category_num;$n1++) {
    $total = 0;
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        $total = $total + @{'num_category'.$n1}[$n2];
    }
    if ($total > 0) {
        $m=1;
        for ($n2=$sample_num-1;$n2>=1;$n2--) {
            @{'freq_category'.$n1}[$m] = @{'num_category'.$n1}[$n2]/$total;
            $m++;
        }
    }
    if ($total == 0) {
        for ($n2=1;$n2<=$sample_num-1;$n2++) {
            @{'freq_category'.$n1}[$n2] = 0;
        }
    }
    @{'freq_category'.$n1}[0] = 0;
    @{'freq_category'.$n1}[$sample_num] = 1;
}

for ($n1=1;$n1<=$mutation_category_num;$n1++) {
    printf "@{'freq_category'.$n1}\n";
}

for ($n1=1;$n1<=$total_seq_num;$n1++) {
    $n2 = $n1-1;
    ${'SEQ'.$n1} = @str1[2*$n2+1];
    @{'SEQ'.$n1} = split //, ${'SEQ'.$n1};
}

$n2 = $sample_start-1;
for ($n1=1;$n1<=$sample_num;$n1++) {
    ${'mSEQ'.$n1} = @str3[2*$n2+1];
    @{'mSEQ'.$n1} = split //, ${'mSEQ'.$n1};
    $n2++;
}


open (OUT1, ">>anc_site_probs.txt");


$n1=0;

for ($n=0;$n<=$#SEQ1-1;$n++) {
    $site_pattern = "";
    for ($x=1;$x<=$total_seq_num;$x++) {
        $site_pattern = $site_pattern.@{'SEQ'.$x}[$n];
    }
    #printf "$site_pattern\n";
    #$site_pattern = "@SEQ1[$n]@SEQ2[$n]@SEQ3[$n]@SEQ4[$n]@SEQ5[$n]@SEQ6[$n]";
    $nref1=0;
    $nref2=0;
    
    $site = $n;
    $ref1 = @mSEQ1[$site];
    $nref1 = 1;
    for ($n2=2;$n2<=$sample_num;$n2++) {
        
        if (@{'mSEQ'.$n2}[$site] eq $ref1) {
            $nref1++;
        }
        if (@{'mSEQ'.$n2}[$site] ne $ref1) {
            $ref2 = @{'mSEQ'.$n2}[$site];
            $nref2++;
        }
        if ((@{'mSEQ'.$n2}[$site] ne $ref1) && (@{'mSEQ'.$n2}[$site] ne $ref2)) {
            printf "error $ref1 $ref2 @{'mSEQ'.$n2}[$site]\n";
        }
    }
    $freq_ref1 = $nref1;
    $freq_ref2 = $nref2;
    if ($nref2 == 0) {
        $ref2 = 'nan';
    }
    $n1++;
    
    #printf "$ref1\t$ref2\t$nref1\t$nref2\n";
    
    foreach $_ (@str2) {
        if ($_ =~ /(\d+)\s+$site_pattern:\s(\w+)/) {
            
            print (OUT1 "$n\t");
            
            $num = $1;
            # find the inferred state at the target node
            $preref = $2;
            @preref = split //, $preref;
            $ref_m = @preref[$target_node-1];
            $ref_ms = @preref[0];
            
            if ($ref_m eq $ref1) {
                $freq = $freq_ref1;
            }
            if ($ref_m eq $ref2) {
                $freq = $freq_ref2;
            }
            
            $total_dp = 0;
            $total = 0;
            # weight the probability using the expected SFS (calculate the total of the weighted probabilities)
            while ($_ =~ /(\w+)\s+\((\d\.\d+)\)/g) {
                $preref = $1;
                
                
                if ($preref ne "total") {
                    @preref = split //, $preref;
                    $ref_m = @preref[$target_node-1];
                    
                    if ($ref_m eq $ref1) {
                        if ($ref_m eq "A") {
                            if ($ref2 eq "nan") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_AT}[$freq_ref1]);
                            }
                            if ($ref2 eq "T") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_AT}[$freq_ref1]);
                            }
                            if ($ref2 eq "G") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_AG}[$freq_ref1]);
                            }
                            if ($ref2 eq "C") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_AC}[$freq_ref1]);
                            }
                        }
                        if ($ref_m eq "T") {
                            if ($ref2 eq "A") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_TA}[$freq_ref1]);
                            }
                            if ($ref2 eq "nan") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_TA}[$freq_ref1]);
                            }
                            if ($ref2 eq "G") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_TG}[$freq_ref1]);
                            }
                            if ($ref2 eq "C") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_TC}[$freq_ref1]);
                            }
                        }
                        if ($ref_m eq "G") {
                            if ($ref2 eq "A") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_GA}[$freq_ref1]);
                            }
                            if ($ref2 eq "T") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_GT}[$freq_ref1]);
                            }
                            if ($ref2 eq "nan") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_GC}[$freq_ref1]);
                            }
                            if ($ref2 eq "C") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_GC}[$freq_ref1]);
                            }
                        }
                        if ($ref_m eq "C") {
                            if ($ref2 eq "A") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_CA}[$freq_ref1]);
                            }
                            if ($ref2 eq "T") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_CT}[$freq_ref1]);
                            }
                            if ($ref2 eq "G") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_CG}[$freq_ref1]);
                            }
                            if ($ref2 eq "nan") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_CG}[$freq_ref1]);
                            }
                        }
                    }
                    
                    if ($ref_m eq $ref2) {
                        if ($ref_m eq "A") {
                            if ($ref1 eq "A") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_AT}[$freq_ref2]);
                            }
                            if ($ref1 eq "T") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_AT}[$freq_ref2]);
                            }
                            if ($ref1 eq "G") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_AG}[$freq_ref2]);
                            }
                            if ($ref1 eq "C") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_AC}[$freq_ref2]);
                            }
                        }
                        if ($ref_m eq "T") {
                            if ($ref1 eq "A") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_TA}[$freq_ref2]);
                            }
                            if ($ref1 eq "T") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_TA}[$freq_ref2]);
                            }
                            if ($ref1 eq "G") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_TG}[$freq_ref2]);
                            }
                            if ($ref1 eq "C") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_TC}[$freq_ref2]);
                            }
                        }
                        if ($ref_m eq "G") {
                            if ($ref1 eq "A") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_GA}[$freq_ref2]);
                            }
                            if ($ref1 eq "T") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_GT}[$freq_ref2]);
                            }
                            if ($ref1 eq "G") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_GC}[$freq_ref2]);
                            }
                            if ($ref1 eq "C") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_GC}[$freq_ref2]);
                            }
                        }
                        if ($ref_m eq "C") {
                            if ($ref1 eq "A") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_CA}[$freq_ref2]);
                            }
                            if ($ref1 eq "T") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_CT}[$freq_ref2]);
                            }
                            if ($ref1 eq "G") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_CG}[$freq_ref2]);
                            }
                            if ($ref1 eq "C") {
                                $total = $total + ($2*@{'freq_category'.$cate_check_CG}[$freq_ref2]);
                            }
                        }
                    }
                    # put probability 0 if the inferred state is different from both of the polymorphic states
                    if (($ref_m ne $ref1) && ($ref_m ne $ref2)) {
                        $total = $total + 0.00;
                        $total_dp = $total_dp + $5;
                    }
                }
            }
            
            
            # weight the probability using the expected SFS (calculate the weighted probability of each ancestral state)
            if ($total_dp <= 1.00) {
                #if ($total == 0) {
                #printf "$site\t$ref_m\t$ref1\t$ref2\n";
                #}
                #print (OUT4 "$site_pattern\n");
                while ($_ =~ /(\w+)\s+\((\d\.\d+)\)/g) {
                    $preref = $1;
                    @preref = split //, $preref;
                    $ref_m = @preref[$target_node-1];
                    
                    if ($ref_m eq $ref1) {
                        if ($ref_m eq "A") {
                            if ($ref2 eq "nan") {
                                $X = ($2*@{'freq_category'.$cate_check_AT}[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "T") {
                                $X = ($2*@{'freq_category'.$cate_check_AT}[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "G") {
                                $X = ($2*@{'freq_category'.$cate_check_AG}[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "C") {
                                $X = ($2*@{'freq_category'.$cate_check_AC}[$freq_ref1])/$total;
                            }
                        }
                        if ($ref_m eq "T") {
                            if ($ref2 eq "A") {
                                $X = ($2*@{'freq_category'.$cate_check_TA}[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "nan") {
                                $X = ($2*@{'freq_category'.$cate_check_TA}[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "G") {
                                $X = ($2*@{'freq_category'.$cate_check_TG}[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "C") {
                                $X = ($2*@{'freq_category'.$cate_check_TC}[$freq_ref1])/$total;
                            }
                        }
                        if ($ref_m eq "G") {
                            if ($ref2 eq "A") {
                                $X = ($2*@{'freq_category'.$cate_check_GA}[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "T") {
                                $X = ($2*@{'freq_category'.$cate_check_GT}[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "nan") {
                                $X = ($2*@{'freq_category'.$cate_check_GC}[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "C") {
                                $X = ($2*@{'freq_category'.$cate_check_GC}[$freq_ref1])/$total;
                            }
                        }
                        if ($ref_m eq "C") {
                            if ($ref2 eq "A") {
                                $X = ($2*@{'freq_category'.$cate_check_CA}[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "T") {
                                $X = ($2*@{'freq_category'.$cate_check_CT}[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "G") {
                                $X = ($2*@{'freq_category'.$cate_check_CG}[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "nan") {
                                $X = ($2*@{'freq_category'.$cate_check_CG}[$freq_ref1])/$total;
                            }
                        }
                    }
                    
                    if ($ref_m eq $ref2) {
                        if ($ref_m eq "A") {
                            if ($ref1 eq "A") {
                                $X = ($2*@{'freq_category'.$cate_check_AT}[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "T") {
                                $X = ($2*@{'freq_category'.$cate_check_AT}[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "G") {
                                $X = ($2*@{'freq_category'.$cate_check_AG}[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "C") {
                                $X = ($2*@{'freq_category'.$cate_check_AC}[$freq_ref2])/$total;
                            }
                        }
                        if ($ref_m eq "T") {
                            if ($ref1 eq "A") {
                                $X = ($2*@{'freq_category'.$cate_check_TA}[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "T") {
                                $X = ($2*@{'freq_category'.$cate_check_TA}[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "G") {
                                $X = ($2*@{'freq_category'.$cate_check_TG}[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "C") {
                                $X = ($2*@{'freq_category'.$cate_check_TC}[$freq_ref2])/$total;
                            }
                        }
                        if ($ref_m eq "G") {
                            if ($ref1 eq "A") {
                                $X = ($2*@{'freq_category'.$cate_check_GA}[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "T") {
                                $X = ($2*@{'freq_category'.$cate_check_GT}[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "G") {
                                $X = ($2*@{'freq_category'.$cate_check_GC}[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "C") {
                                $X = ($2*@{'freq_category'.$cate_check_GC}[$freq_ref2])/$total;
                            }
                        }
                        if ($ref_m eq "C") {
                            if ($ref1 eq "A") {
                                $X = ($2*@{'freq_category'.$cate_check_CA}[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "T") {
                                $X = ($2*@{'freq_category'.$cate_check_CT}[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "G") {
                                $X = ($2*@{'freq_category'.$cate_check_CG}[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "C") {
                                $X = ($2*@{'freq_category'.$cate_check_CG}[$freq_ref2])/$total;
                            }
                        }
                    }
                    
                    if (($ref_m ne $ref1) && ($ref_m ne $ref2)) {
                        $X = 0.00;
                    }
                    
                    print (OUT1 "$ref_m\t$X\t");
                    
                    if ($ref_ms ne $ref_m) {
                        if ($ref_ms eq "T") {
                            if ($ref_m eq "C") {$num_TC = $num_TC + $X;}
                            if ($ref_m eq "A") {$num_TA = $num_TA + $X;}
                            if ($ref_m eq "G") {$num_TG = $num_TG + $X;}
                        }
                        if ($ref_ms eq "C") {
                            if ($ref_m eq "T") {$num_CT = $num_CT + $X;}
                            if ($ref_m eq "A") {$num_CA = $num_CA + $X;}
                            if ($ref_m eq "G") {$num_CG = $num_CG + $X;}
                        }
                        if ($ref_ms eq "A") {
                            if ($ref_m eq "C") {$num_AC = $num_AC + $X;}
                            if ($ref_m eq "T") {$num_AT = $num_AT + $X;}
                            if ($ref_m eq "G") {$num_AG = $num_AG + $X;}
                        }
                        if ($ref_ms eq "G") {
                            if ($ref_m eq "C") {$num_GC = $num_GC + $X;}
                            if ($ref_m eq "A") {$num_GA = $num_GA + $X;}
                            if ($ref_m eq "T") {$num_GT = $num_GT + $X;}
                        }
                    }
                }
            }
            if ($total_dp > 1.00) {
                #print (OUT5 "$site_pattern\n");
                #print (OUT6 "$n\n");
                while ($_ =~ /(\w)(\w)(\w)(\w)\s+\((\d\.\d+)\)/g) {
                    
                    print (OUT1 "$ref_m\tnan\t");
                    
                }
            }
            
            print (OUT1 "\n");
        }
    }
}


open (OUT7, ">>AWP_substitution_for_post_weighted_collapse.txt");
print (OUT7 "TC\tTA\tTG\tCT\tCA\tCG\tAT\tAC\tAG\tGT\tGC\tGA\n");
print (OUT7 "$num_TC\t$num_TA\t$num_TG\t$num_CT\t$num_CA\t$num_CG\t$num_AT\t$num_AC\t$num_AG\t$num_GT\t$num_GC\t$num_GA\n");
close (OUT7);
