open (IN1,"sample_seq_m.txt");
@str1 = <IN1>;
open (IN2,"BASEML_input_melpoly.MFA");
@str2 = <IN2>;
open (IN3,"site_used_for_collapse_method.txt");
@str3 = <IN3>;
open (IN4,"site_not_used_for_higher_than_0%_ms_analysis.txt");
@str4 = <IN4>;


$n=1;
foreach $_ (@str1) {
    if ($_ =~ /^\w\w\w\w\w\w\w/) {
        ${'SEQ'.$n} = $_;
        @{'SEQ'.$n} = split //, ${'SEQ'.$n};
        $n++;
    }
}

$ancseq_m = @str2[1];
$ancseq_t = @str2[5];
$ancseq_e = @str2[9];

@ANCSEQ_m = split //, $ancseq_m;
@ANCSEQ_t = split //, $ancseq_t;
@ANCSEQ_e = split //, $ancseq_e;

for ($n1=1;$n1<=10;$n1++) {
    @freq_spec_TC[$n1]=0;
    @freq_spec_TA[$n1]=0;
    @freq_spec_TG[$n1]=0;
    @freq_spec_CT[$n1]=0;
    @freq_spec_CA[$n1]=0;
    @freq_spec_CG[$n1]=0;
    @freq_spec_AT[$n1]=0;
    @freq_spec_AC[$n1]=0;
    @freq_spec_AG[$n1]=0;
    @freq_spec_GT[$n1]=0;
    @freq_spec_GC[$n1]=0;
    @freq_spec_GA[$n1]=0;
}

$n=0;
foreach $_ (@str3) {
    $AAAA=0;
    $site = $_;
    $der1 = "X";
    $der2 = "Y";
    
    for ($n1=1;$n1<=10;$n1++) {
        ${'sample'.$n1} = @{'SEQ'.$n1}[$site];
    }
    
    #use ancestral sequence after the filtering of >2 states sites
    $anc_site = @ANCSEQ_m[$n];
    foreach $_ (@str4) {
        if ($n == $_) {
            $AAAA=1;
        }
    }
    $n++;
    
    if ($AAAA==0) {
        $check1=0;
        $check2=0;
        
        for ($n1=1;$n1<=10;$n1++) {
            if (${'sample'.$n1} ne $anc_site) {
                if (($der1 eq "X") || (${'sample'.$n1} eq $der1)) {
                    $check1++;
                    $der1 = ${'sample'.$n1};
                }
                elsif (${'sample'.$n1} ne $der1) {
                    $check2++;
                    $der2 = ${'sample'.$n1};
                }
            }
        }
        
        if ($anc_site eq "T") {
            if ($der1 eq "C") {
                @freq_spec_TC[$check1] = @freq_spec_TC[$check1] + 1;
            }
            if ($der1 eq "A") {
                @freq_spec_TA[$check1] = @freq_spec_TA[$check1] + 1;
            }
            if ($der1 eq "G") {
                @freq_spec_TG[$check1] = @freq_spec_TG[$check1] + 1;
            }
            if ($der2 eq "C") {
                @freq_spec_TC[$check2] = @freq_spec_TC[$check2] + 1;
            }
            if ($der2 eq "A") {
                @freq_spec_TA[$check2] = @freq_spec_TA[$check2] + 1;
            }
            if ($der2 eq "G") {
                @freq_spec_TG[$check2] = @freq_spec_TG[$check2] + 1;
            }
        }
        if ($anc_site eq "C") {
            if ($der1 eq "T") {
                @freq_spec_CT[$check1] = @freq_spec_CT[$check1] + 1;
            }
            if ($der1 eq "A") {
                @freq_spec_CA[$check1] = @freq_spec_CA[$check1] + 1;
            }
            if ($der1 eq "G") {
                @freq_spec_CG[$check1] = @freq_spec_CG[$check1] + 1;
            }
            if ($der2 eq "T") {
                @freq_spec_CT[$check2] = @freq_spec_CT[$check2] + 1;
            }
            if ($der2 eq "A") {
                @freq_spec_CA[$check2] = @freq_spec_CA[$check2] + 1;
            }
            if ($der2 eq "G") {
                @freq_spec_CG[$check2] = @freq_spec_CG[$check2] + 1;
            }
        }
        if ($anc_site eq "A") {
            if ($der1 eq "T") {
                @freq_spec_AT[$check1] = @freq_spec_AT[$check1] + 1;
            }
            if ($der1 eq "C") {
                @freq_spec_AC[$check1] = @freq_spec_AC[$check1] + 1;
            }
            if ($der1 eq "G") {
                @freq_spec_AG[$check1] = @freq_spec_AG[$check1] + 1;
            }
            if ($der2 eq "T") {
                @freq_spec_AT[$check2] = @freq_spec_AT[$check2] + 1;
            }
            if ($der2 eq "C") {
                @freq_spec_AC[$check2] = @freq_spec_AC[$check2] + 1;
            }
            if ($der2 eq "G") {
                @freq_spec_AG[$check2] = @freq_spec_AG[$check2] + 1;
            }
        }
        if ($anc_site eq "G") {
            if ($der1 eq "T") {
                @freq_spec_GT[$check1] = @freq_spec_GT[$check1] + 1;
            }
            if ($der1 eq "C") {
                @freq_spec_GC[$check1] = @freq_spec_GC[$check1] + 1;
            }
            if ($der1 eq "A") {
                @freq_spec_GA[$check1] = @freq_spec_GA[$check1] + 1;
            }
            if ($der2 eq "T") {
                @freq_spec_GT[$check2] = @freq_spec_GT[$check2] + 1;
            }
            if ($der2 eq "C") {
                @freq_spec_GC[$check2] = @freq_spec_GC[$check2] + 1;
            }
            if ($der2 eq "A") {
                @freq_spec_GA[$check2] = @freq_spec_GA[$check2] + 1;
            }
        }
    }
}

open (OUT1, ">>actual_frequency_spectrum_m_0%_ms.txt");
print (OUT1 "TC\tTA\tTG\tCT\tCA\tCG\tAT\tAC\tAG\tGT\tGC\tGA\n");
close (OUT1);
for ($n1=1;$n1<=10;$n1++) {
    open (OUT1, ">>actual_frequency_spectrum_m_0%_ms.txt");
    print (OUT1 "@freq_spec_TC[$n1]\t@freq_spec_TA[$n1]\t@freq_spec_TG[$n1]\t@freq_spec_CT[$n1]\t@freq_spec_CA[$n1]\t@freq_spec_CG[$n1]\t@freq_spec_AT[$n1]\t@freq_spec_AC[$n1]\t@freq_spec_AG[$n1]\t@freq_spec_GT[$n1]\t@freq_spec_GC[$n1]\t@freq_spec_GA[$n1]\n");
    close (OUT1);
}


