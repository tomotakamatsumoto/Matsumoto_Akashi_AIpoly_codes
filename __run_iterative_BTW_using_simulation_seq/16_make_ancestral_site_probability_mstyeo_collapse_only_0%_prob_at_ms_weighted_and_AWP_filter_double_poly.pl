open (IN1,"collapse.MFA.mfa");
@str1 = <IN1>;
open (IN2,"fltrst");
@str2 = <IN2>;
open (IN3,"sample_seq_m.txt");
@str3 = <IN3>;
open (IN4,"site_used_for_collapse_method.txt");
@str4 = <IN4>;

#@nPRF = (0.000, 0.040, 0.045, 0.051, 0.059, 0.071, 0.089, 0.119, 0.178, 0.348, 1.000);
#@nPRF = (0.000, 0.040, 0.045, 0.051, 0.059, 0.071, 0.089, 0.119, 0.178, 0.348, 1.000);
#@nPRF = (0.000, 0.040, 0.045, 0.051, 0.059, 0.071, 0.089, 0.119, 0.178, 0.348, 1.000);

@nPRF = (0.000, 0.040, 0.045, 0.051, 0.059, 0.071, 0.089, 0.119, 0.178, 0.348, 1.000);
@upPRF = (0.000, 0.041, 0.046, 0.052, 0.060, 0.072, 0.089, 0.118, 0.176, 0.347, 1.000);
@puPRF = (0.000, 0.038, 0.043, 0.049, 0.058, 0.070, 0.088, 0.118, 0.179, 0.358, 1.000);


$SEQ1 = @str1[1];
$SEQ2 = @str1[3];
$SEQ3 = @str1[5];
$SEQ4 = @str1[7];
$SEQ5 = @str1[9];
$SEQ6 = @str1[11];
$SEQ7 = @str1[13];
$SEQ8 = @str1[15];
$SEQ9 = @str1[17];
$SEQ10 = @str1[19];
$SEQ11 = @str1[21];
$SEQ12 = @str1[23];

@SEQ1 = split //, $SEQ1;
@SEQ2 = split //, $SEQ2;
@SEQ3 = split //, $SEQ3;
@SEQ4 = split //, $SEQ4;
@SEQ5 = split //, $SEQ5;
@SEQ6 = split //, $SEQ6;
@SEQ7 = split //, $SEQ7;
@SEQ8 = split //, $SEQ8;
@SEQ9 = split //, $SEQ9;
@SEQ10 = split //, $SEQ10;
@SEQ11 = split //, $SEQ11;
@SEQ12 = split //, $SEQ12;

$mSEQ1 = @str3[1];
$mSEQ2 = @str3[3];
$mSEQ3 = @str3[5];
$mSEQ4 = @str3[7];
$mSEQ5 = @str3[9];
$mSEQ6 = @str3[11];
$mSEQ7 = @str3[13];
$mSEQ8 = @str3[15];
$mSEQ9 = @str3[17];
$mSEQ10 = @str3[19];

@mSEQ1 = split //, $mSEQ1;
@mSEQ2 = split //, $mSEQ2;
@mSEQ3 = split //, $mSEQ3;
@mSEQ4 = split //, $mSEQ4;
@mSEQ5 = split //, $mSEQ5;
@mSEQ6 = split //, $mSEQ6;
@mSEQ7 = split //, $mSEQ7;
@mSEQ8 = split //, $mSEQ8;
@mSEQ9 = split //, $mSEQ9;
@mSEQ10 = split //, $mSEQ10;


open (OUT1, ">>anc_site_prob_m_0%_ms.txt");
#open (OUT2, ">>anc_site_prob_t.txt");
#open (OUT3, ">>anc_site_prob_e.txt");
open (OUT4, ">>site_pattern_with_higher_than_0%_ms_node.txt");
open (OUT5, ">>site_pattern_with_lower_than_0%_ms_node.txt");
open (OUT6, ">>site_not_used_for_higher_than_0%_ms_analysis.txt");

$n1=0;

for ($n=0;$n<=$#SEQ1-1;$n++) {
    $site_pattern = "@SEQ1[$n]@SEQ2[$n]@SEQ3[$n]@SEQ4[$n]@SEQ5[$n]@SEQ6[$n]@SEQ7[$n]@SEQ8[$n]@SEQ9[$n]@SEQ10[$n]@SEQ11[$n]@SEQ12[$n]";
    $nref1=0;
    $nref2=0;
    
    
    $site = @str4[$n1];
    #printf "@SEQ1[$n]@SEQ2[$n] @mSEQ1[$site]@mSEQ2[$site]@mSEQ3[$site]@mSEQ4[$site]@mSEQ5[$site]@mSEQ6[$site]@mSEQ7[$site]@mSEQ8[$site]@mSEQ9[$site]@mSEQ10[$site]\n";
    $ref1 = @mSEQ1[$site];
    $nref1 = 1;
    for ($n2=2;$n2<=10;$n2++) {
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
        if ($_ =~ /(\d+)\s+$site_pattern:\s\w(\w)/) {
            
            print (OUT1 "$n\t");
            #print (OUT2 "$n\t");
            #print (OUT3 "$n\t");
            $num = $1;
            $ref_m = $2;
            if ($ref_m eq $ref1) {
                $freq = $freq_ref1;
            }
            if ($ref_m eq $ref2) {
                $freq = $freq_ref2;
            }
            
            $total_dp = 0;
            $total = 0;
            while ($_ =~ /(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\w)\s+\((\d\.\d+)\)/g) {
                if ($3 eq $ref1) {
                    if ($3 eq "A") {
                        if ($ref2 eq "nan") {
                            $total = $total + ($11*@nPRF[$freq_ref1]);
                        }
                        if ($ref2 eq "T") {
                            $total = $total + ($11*@nPRF[$freq_ref1]);
                        }
                        if ($ref2 eq "G") {
                            $total = $total + ($11*@upPRF[$freq_ref1]);
                        }
                        if ($ref2 eq "C") {
                            $total = $total + ($11*@upPRF[$freq_ref1]);
                        }
                    }
                    if ($3 eq "T") {
                        if ($ref2 eq "A") {
                            $total = $total + ($11*@nPRF[$freq_ref1]);
                        }
                        if ($ref2 eq "nan") {
                            $total = $total + ($11*@nPRF[$freq_ref1]);
                        }
                        if ($ref2 eq "G") {
                            $total = $total + ($11*@upPRF[$freq_ref1]);
                        }
                        if ($ref2 eq "C") {
                            $total = $total + ($11*@upPRF[$freq_ref1]);
                        }
                    }
                    if ($3 eq "G") {
                        if ($ref2 eq "A") {
                            $total = $total + ($11*@puPRF[$freq_ref1]);
                        }
                        if ($ref2 eq "T") {
                            $total = $total + ($11*@puPRF[$freq_ref1]);
                        }
                        if ($ref2 eq "nan") {
                            $total = $total + ($11*@nPRF[$freq_ref1]);
                        }
                        if ($ref2 eq "C") {
                            $total = $total + ($11*@nPRF[$freq_ref1]);
                        }
                    }
                    if ($3 eq "C") {
                        if ($ref2 eq "A") {
                            $total = $total + ($11*@puPRF[$freq_ref1]);
                        }
                        if ($ref2 eq "T") {
                            $total = $total + ($11*@puPRF[$freq_ref1]);
                        }
                        if ($ref2 eq "G") {
                            $total = $total + ($11*@nPRF[$freq_ref1]);
                        }
                        if ($ref2 eq "nan") {
                            $total = $total + ($11*@nPRF[$freq_ref1]);
                        }
                    }
                }
                
                if ($3 eq $ref2) {
                    if ($3 eq "A") {
                        if ($ref1 eq "A") {
                            $total = $total + ($11*@nPRF[$freq_ref2]);
                        }
                        if ($ref1 eq "T") {
                            $total = $total + ($11*@nPRF[$freq_ref2]);
                        }
                        if ($ref1 eq "G") {
                            $total = $total + ($11*@upPRF[$freq_ref2]);
                        }
                        if ($ref1 eq "C") {
                            $total = $total + ($11*@upPRF[$freq_ref2]);
                        }
                    }
                    if ($3 eq "T") {
                        if ($ref1 eq "A") {
                            $total = $total + ($11*@nPRF[$freq_ref2]);
                        }
                        if ($ref1 eq "T") {
                            $total = $total + ($11*@nPRF[$freq_ref2]);
                        }
                        if ($ref1 eq "G") {
                            $total = $total + ($11*@upPRF[$freq_ref2]);
                        }
                        if ($ref1 eq "C") {
                            $total = $total + ($11*@upPRF[$freq_ref2]);
                        }
                    }
                    if ($3 eq "G") {
                        if ($ref1 eq "A") {
                            $total = $total + ($11*@puPRF[$freq_ref2]);
                        }
                        if ($ref1 eq "T") {
                            $total = $total + ($11*@puPRF[$freq_ref2]);
                        }
                        if ($ref1 eq "G") {
                            $total = $total + ($11*@nPRF[$freq_ref2]);
                        }
                        if ($ref1 eq "C") {
                            $total = $total + ($11*@nPRF[$freq_ref2]);
                        }
                    }
                    if ($3 eq "C") {
                        if ($ref1 eq "A") {
                            $total = $total + ($11*@puPRF[$freq_ref2]);
                        }
                        if ($ref1 eq "T") {
                            $total = $total + ($11*@puPRF[$freq_ref2]);
                        }
                        if ($ref1 eq "G") {
                            $total = $total + ($11*@nPRF[$freq_ref2]);
                        }
                        if ($ref1 eq "C") {
                            $total = $total + ($11*@nPRF[$freq_ref2]);
                        }
                    }
                }
                
                if (($3 ne $ref1) && ($3 ne $ref2)) {
                    $total = $total + 0.00;
                    $total_dp = $total_dp + $11;
                }
            }
            
            
            if ($total_dp <= 1.00) {
                print (OUT4 "$site_pattern\n");
                while ($_ =~ /(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\w)\s+\((\d\.\d+)\)/g) {
                    
                    if ($3 eq $ref1) {
                        if ($3 eq "A") {
                            if ($ref2 eq "nan") {
                                $X = ($11*@nPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "T") {
                                $X = ($11*@nPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "G") {
                                $X = ($11*@upPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "C") {
                                $X = ($11*@upPRF[$freq_ref1])/$total;
                            }
                        }
                        if ($3 eq "T") {
                            if ($ref2 eq "A") {
                                $X = ($11*@nPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "nan") {
                                $X = ($11*@nPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "G") {
                                $X = ($11*@upPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "C") {
                                $X = ($11*@upPRF[$freq_ref1])/$total;
                            }
                        }
                        if ($3 eq "G") {
                            if ($ref2 eq "A") {
                                $X = ($11*@puPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "T") {
                                $X = ($11*@puPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "nan") {
                                $X = ($11*@nPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "C") {
                                $X = ($11*@nPRF[$freq_ref1])/$total;
                            }
                        }
                        if ($3 eq "C") {
                            if ($ref2 eq "A") {
                                $X = ($11*@puPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "T") {
                                $X = ($11*@puPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "G") {
                                $X = ($11*@nPRF[$freq_ref1])/$total;
                            }
                            if ($ref2 eq "nan") {
                                $X = ($11*@nPRF[$freq_ref1])/$total;
                            }
                        }
                    }
                    
                    if ($3 eq $ref2) {
                        if ($3 eq "A") {
                            if ($ref1 eq "A") {
                                $X = ($11*@nPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "T") {
                                $X = ($11*@nPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "G") {
                                $X = ($11*@upPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "C") {
                                $X = ($11*@upPRF[$freq_ref2])/$total;
                            }
                        }
                        if ($3 eq "T") {
                            if ($ref1 eq "A") {
                                $X = ($11*@nPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "T") {
                                $X = ($11*@nPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "G") {
                                $X = ($11*@upPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "C") {
                                $X = ($11*@upPRF[$freq_ref2])/$total;
                            }
                        }
                        if ($3 eq "G") {
                            if ($ref1 eq "A") {
                                $X = ($11*@puPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "T") {
                                $X = ($11*@puPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "G") {
                                $X = ($11*@nPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "C") {
                                $X = ($11*@nPRF[$freq_ref2])/$total;
                            }
                        }
                        if ($3 eq "C") {
                            if ($ref1 eq "A") {
                                $X = ($11*@puPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "T") {
                                $X = ($11*@puPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "G") {
                                $X = ($11*@nPRF[$freq_ref2])/$total;
                            }
                            if ($ref1 eq "C") {
                                $X = ($11*@nPRF[$freq_ref2])/$total;
                            }
                        }
                    }
                    
                    if (($3 ne $ref1) && ($3 ne $ref2)) {
                        $X = 0.00;
                    }
                    
                    print (OUT1 "$3\t$X\t");
                    #print (OUT2 "$6\t$X\t");
                    #print (OUT3 "$9\t$X\t");
                    
                    if ($2 ne $3) {
                        if ($2 eq "T") {
                            if ($3 eq "C") {$num_TC = $num_TC + $X;}
                            if ($3 eq "A") {$num_TA = $num_TA + $X;}
                            if ($3 eq "G") {$num_TG = $num_TG + $X;}
                        }
                        if ($2 eq "C") {
                            if ($3 eq "T") {$num_CT = $num_CT + $X;}
                            if ($3 eq "A") {$num_CA = $num_CA + $X;}
                            if ($3 eq "G") {$num_CG = $num_CG + $X;}
                        }
                        if ($2 eq "A") {
                            if ($3 eq "C") {$num_AC = $num_AC + $X;}
                            if ($3 eq "T") {$num_AT = $num_AT + $X;}
                            if ($3 eq "G") {$num_AG = $num_AG + $X;}
                        }
                        if ($2 eq "G") {
                            if ($3 eq "C") {$num_GC = $num_GC + $X;}
                            if ($3 eq "A") {$num_GA = $num_GA + $X;}
                            if ($3 eq "T") {$num_GT = $num_GT + $X;}
                        }
                    }
                }
            }
            if ($total_dp > 1.00) {
                print (OUT5 "$site_pattern\n");
                print (OUT6 "$n\n");
                while ($_ =~ /(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\w)\s+\((\d\.\d+)\)/g) {
                    
                    print (OUT1 "$3\tnan\t");
                    #print (OUT2 "$6\tnan\t");
                    #print (OUT3 "$9\tnan\t");
                }
            }
            
            print (OUT1 "\n");
            #print (OUT2 "\n");
            #print (OUT3 "\n");
        }
    }
}
open (OUT7, ">>AWP_substitution_for_post_weighted_collapse.txt");
print (OUT7 "TC\tTA\tTG\tCT\tCA\tCG\tAT\tAC\tAG\tGT\tGC\tGA\n");
print (OUT7 "$num_TC\t$num_TA\t$num_TG\t$num_CT\t$num_CA\t$num_CG\t$num_AT\t$num_AC\t$num_AG\t$num_GT\t$num_GC\t$num_GA\n");
close (OUT7);
