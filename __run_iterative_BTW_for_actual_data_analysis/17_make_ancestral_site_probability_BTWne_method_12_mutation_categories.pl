open (IN1,"collapse.MFA.mfa");
@str1 = <IN1>;
open (IN2,"fltrst");
@str2 = <IN2>;
open (IN3,"sample_seq.txt");
@str3 = <IN3>;

# read integers : total number of sequences (collapse seqs + outgroups), number of alleles in the focused species, number of internal nodes in the phylogeny, order of the target node to which the ancestral state is inferred
$total_seq_num = $ARGV[0];
$num_internal_nodes = $ARGV[1];
$target_node = $ARGV[2];
$sample_start = $ARGV[3];
$sample_end = $ARGV[4];

$sample_num = $sample_end - $sample_start + 1;

# make SFSne
$sum = 0;
for ($n=1;$n<=$sample_num-1;$n++) {
    $sum = $sum + (1/$n);
}
$m=1;
for ($n=$sample_num-1;$n>=1;$n--) {
    @TCPRF[$m] = (1/$n)/$sum;
    @TAPRF[$m] = (1/$n)/$sum;
    @TGPRF[$m] = (1/$n)/$sum;
    
    @CTPRF[$m] = (1/$n)/$sum;
    @CAPRF[$m] = (1/$n)/$sum;
    @CGPRF[$m] = (1/$n)/$sum;
    
    @ATPRF[$m] = (1/$n)/$sum;
    @ACPRF[$m] = (1/$n)/$sum;
    @AGPRF[$m] = (1/$n)/$sum;
    
    @GTPRF[$m] = (1/$n)/$sum;
    @GCPRF[$m] = (1/$n)/$sum;
    @GAPRF[$m] = (1/$n)/$sum;
    $m++;
}

@TCPRF[$sample_num] = 1;
@TAPRF[$sample_num] = 1;
@TGPRF[$sample_num] = 1;
@CTPRF[$sample_num] = 1;
@CAPRF[$sample_num] = 1;
@CGPRF[$sample_num] = 1;
@ATPRF[$sample_num] = 1;
@ACPRF[$sample_num] = 1;
@AGPRF[$sample_num] = 1;
@GTPRF[$sample_num] = 1;
@GCPRF[$sample_num] = 1;
@GAPRF[$sample_num] = 1;

@TCPRF[0] = 0;
@TAPRF[0] = 0;
@TGPRF[0] = 0;
@CTPRF[0] = 0;
@CAPRF[0] = 0;
@CGPRF[0] = 0;
@ATPRF[0] = 0;
@ACPRF[0] = 0;
@AGPRF[0] = 0;
@GTPRF[0] = 0;
@GCPRF[0] = 0;
@GAPRF[0] = 0;


for ($n1=1;$n1<=$total_seq_num;$n1++) {
    $n2 = $n1-1;
    ${'SEQ'.$n1} = @str1[2*$n2+1];
    @{'SEQ'.$n1} = split //, ${'SEQ'.$n1};
}

for ($n1=1;$n1<=$sample_num;$n1++) {
    $n2 = $sample_start-1;
    ${'mSEQ'.$n1} = @str3[2*$n2+1];
    @{'mSEQ'.$n1} = split //, ${'mSEQ'.$n1};
}


open (OUT1, ">>anc_site_probs.txt");
#open (OUT4, ">>site_pattern_with_higher_than_0%_ms_node.txt");
#open (OUT5, ">>site_pattern_with_lower_than_0%_ms_node.txt");
#open (OUT6, ">>site_not_used_for_higher_than_0%_ms_analysis.txt");

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
                                $total = $total + ($2*@ATPRF[$freq_ref1]);
                            }
                            if ($ref2 eq "T") {
                                $total = $total + ($2*@ATPRF[$freq_ref1]);
                            }
                            if ($ref2 eq "G") {
                                $total = $total + ($2*@AGPRF[$freq_ref1]);
                            }
                            if ($ref2 eq "C") {
                                $total = $total + ($2*@ACPRF[$freq_ref1]);
                            }
                        }
                        if ($ref_m eq "T") {
                            if ($ref2 eq "A") {
                                $total = $total + ($2*@TAPRF[$freq_ref1]);
                            }
                            if ($ref2 eq "nan") {
                                $total = $total + ($2*@TAPRF[$freq_ref1]);
                            }
                            if ($ref2 eq "G") {
                                $total = $total + ($2*@TGPRF[$freq_ref1]);
                            }
                            if ($ref2 eq "C") {
                                $total = $total + ($2*@TCPRF[$freq_ref1]);
                            }
                        }
                        if ($ref_m eq "G") {
                            if ($ref2 eq "A") {
                                $total = $total + ($2*@GAPRF[$freq_ref1]);
                            }
                            if ($ref2 eq "T") {
                                $total = $total + ($2*@GTPRF[$freq_ref1]);
                            }
                            if ($ref2 eq "nan") {
                                $total = $total + ($2*@GCPRF[$freq_ref1]);
                            }
                            if ($ref2 eq "C") {
                                $total = $total + ($2*@GCPRF[$freq_ref1]);
                            }
                        }
                        if ($ref_m eq "C") {
                            if ($ref2 eq "A") {
                                $total = $total + ($2*@CAPRF[$freq_ref1]);
                            }
                            if ($ref2 eq "T") {
                                $total = $total + ($2*@CTPRF[$freq_ref1]);
                            }
                            if ($ref2 eq "G") {
                                $total = $total + ($2*@CGPRF[$freq_ref1]);
                            }
                            if ($ref2 eq "nan") {
                                $total = $total + ($2*@CGPRF[$freq_ref1]);
                            }
                        }
                    }
                    
                    if ($ref_m eq $ref2) {
                        if ($ref_m eq "A") {
                            if ($ref1 eq "A") {
                                $total = $total + ($2*@ATPRF[$freq_ref2]);
                            }
                            if ($ref1 eq "T") {
                                $total = $total + ($2*@ATPRF[$freq_ref2]);
                            }
                            if ($ref1 eq "G") {
                                $total = $total + ($2*@AGPRF[$freq_ref2]);
                            }
                            if ($ref1 eq "C") {
                                $total = $total + ($2*@ACPRF[$freq_ref2]);
                            }
                        }
                        if ($ref_m eq "T") {
                            if ($ref1 eq "A") {
                                $total = $total + ($2*@TAPRF[$freq_ref2]);
                            }
                            if ($ref1 eq "T") {
                                $total = $total + ($2*@TAPRF[$freq_ref2]);
                            }
                            if ($ref1 eq "G") {
                                $total = $total + ($2*@TGPRF[$freq_ref2]);
                            }
                            if ($ref1 eq "C") {
                                $total = $total + ($2*@TCPRF[$freq_ref2]);
                            }
                        }
                        if ($ref_m eq "G") {
                            if ($ref1 eq "A") {
                                $total = $total + ($2*@GAPRF[$freq_ref2]);
                            }
                            if ($ref1 eq "T") {
                                $total = $total + ($2*@GTPRF[$freq_ref2]);
                            }
                            if ($ref1 eq "G") {
                                $total = $total + ($2*@GCPRF[$freq_ref2]);
                            }
                            if ($ref1 eq "C") {
                                $total = $total + ($2*@GCPRF[$freq_ref2]);
                            }
                        }
                        if ($ref_m eq "C") {
                            if ($ref1 eq "A") {
                                $total = $total + ($2*@CAPRF[$freq_ref2]);
                            }
                            if ($ref1 eq "T") {
                                $total = $total + ($2*@CTPRF[$freq_ref2]);
                            }
                            if ($ref1 eq "G") {
                                $total = $total + ($2*@CGPRF[$freq_ref2]);
                            }
                            if ($ref1 eq "C") {
                                $total = $total + ($2*@CGPRF[$freq_ref2]);
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
                #print (OUT4 "$site_pattern\n");
                while ($_ =~ /(\w+)\s+\((\d\.\d+)\)/g) {
                    $preref = $1;
                    @preref = split //, $preref;
                    $ref_m = @preref[$target_node-1];
                    
                    if ($ref_m eq $ref1) {
                        if ($ref_m eq "A") {
                            if ($ref2 eq "nan") {
                                $X = ($2*@ATPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "T") {
                                $X = ($2*@ATPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "G") {
                                $X = ($2*@AGPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "C") {
                                $X = ($2*@ACPRF[$freq_ref1])/$total;
                            }
                        }
                        if ($ref_m eq "T") {
                            if ($ref2 eq "A") {
                                $X = ($2*@TAPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "nan") {
                                $X = ($2*@TAPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "G") {
                                $X = ($2*@TGPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "C") {
                                $X = ($2*@TCPRF[$freq_ref1])/$total;
                            }
                        }
                        if ($ref_m eq "G") {
                            if ($ref2 eq "A") {
                                $X = ($2*@GAPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "T") {
                                $X = ($2*@GTPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "nan") {
                                $X = ($2*@GCPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "C") {
                                $X = ($2*@GCPRF[$freq_ref1])/$total;
                            }
                        }
                        if ($ref_m eq "C") {
                            if ($ref2 eq "A") {
                                $X = ($2*@CAPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "T") {
                                $X = ($2*@CTPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "G") {
                                $X = ($2*@CGPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "nan") {
                                $X = ($2*@CGPRF[$freq_ref1])/$total;
                            }
                        }
                    }
                    
                    if ($ref_m eq $ref2) {
                        if ($ref_m eq "A") {
                            if ($ref1 eq "A") {
                                $X = ($2*@ATPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "T") {
                                $X = ($2*@ATPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "G") {
                                $X = ($2*@AGPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "C") {
                                $X = ($2*@ACPRF[$freq_ref2])/$total;
                            }
                        }
                        if ($ref_m eq "T") {
                            if ($ref1 eq "A") {
                                $X = ($2*@TAPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "T") {
                                $X = ($2*@TAPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "G") {
                                $X = ($2*@TGPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "C") {
                                $X = ($2*@TCPRF[$freq_ref2])/$total;
                            }
                        }
                        if ($ref_m eq "G") {
                            if ($ref1 eq "A") {
                                $X = ($2*@GAPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "T") {
                                $X = ($2*@GTPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "G") {
                                $X = ($2*@GCPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "C") {
                                $X = ($2*@GCPRF[$freq_ref2])/$total;
                            }
                        }
                        if ($ref_m eq "C") {
                            if ($ref1 eq "A") {
                                $X = ($2*@CAPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "T") {
                                $X = ($2*@CTPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "G") {
                                $X = ($2*@CGPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "C") {
                                $X = ($2*@CGPRF[$freq_ref2])/$total;
                            }
                        }
                    }
                    
                    if (($ref_m ne $ref1) && ($ref_m ne $ref2)) {
                        $X = 0.00;
                    }
                    
                    print (OUT1 "$ref_m\t$X\t");
                    
                    #if ($1 ne $ref_m) {
                    #    if ($1 eq "T") {
                    #        if ($3 eq "C") {$num_TC = $num_TC + $X;}
                    #        if ($3 eq "A") {$num_TA = $num_TA + $X;}
                    #        if ($3 eq "G") {$num_TG = $num_TG + $X;}
                    #    }
                    #    if ($1 eq "C") {
                    #        if ($3 eq "T") {$num_CT = $num_CT + $X;}
                    #        if ($3 eq "A") {$num_CA = $num_CA + $X;}
                    #        if ($3 eq "G") {$num_CG = $num_CG + $X;}
                    #    }
                    #    if ($1 eq "A") {
                    #        if ($3 eq "C") {$num_AC = $num_AC + $X;}
                    #        if ($3 eq "T") {$num_AT = $num_AT + $X;}
                    #        if ($3 eq "G") {$num_AG = $num_AG + $X;}
                    #    }
                    #    if ($1 eq "G") {
                    #        if ($3 eq "C") {$num_GC = $num_GC + $X;}
                    #        if ($3 eq "A") {$num_GA = $num_GA + $X;}
                    #        if ($3 eq "T") {$num_GT = $num_GT + $X;}
                    #    }
                    #}
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

#$nnum_TC = $num_TC/($freqT+$freqA);
#$nnum_TA = $num_TA/($freqT+$freqA);
#$nnum_TG = $num_TG/($freqT+$freqA);
#$nnum_CT = $num_CT/($freqC+$freqG);
#$nnum_CA = $num_CA/($freqC+$freqG);
#$nnum_CG = $num_CG/($freqC+$freqG);
#$nnum_AT = $num_AT/($freqT+$freqA);
#$nnum_AC = $num_AC/($freqT+$freqA);
#$nnum_AG = $num_AG/($freqT+$freqA);
#$nnum_GT = $num_GT/($freqC+$freqG);
#$nnum_GC = $num_GC/($freqC+$freqG);
#$nnum_GA = $num_GA/($freqC+$freqG);

#open (OUT7, ">>AWP_substitution_for_post_weighted_collapse.txt");
#print (OUT7 "TC\tTA\tTG\tCT\tCA\tCG\tAT\tAC\tAG\tGT\tGC\tGA\n");
#print (OUT7 "$num_TC\t$num_TA\t$num_TG\t$num_CT\t$num_CA\t$num_CG\t$num_AT\t$num_AC\t$num_AG\t$num_GT\t$num_GC\t$num_GA\n");
#close (OUT7);
