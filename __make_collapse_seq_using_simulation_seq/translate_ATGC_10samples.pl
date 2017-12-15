

#convert current 0~3 nucleotide seq to ATGC seq
open (IN1, "result_m.txt");
@str1 = <IN1>;
$count=1;
$Z=1;
for ($n1=0;$n1<=1000;$n1++) {
    @check[$n1]=1;
}
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /whole\ssequence/) {
        while ($count<=10) {
            $rand = int (rand 1000);
            if (@check[$rand]==1) {
                if (@str1[$n2+$rand+1] =~ /\=\s\"(\d\d\d\d\d+)/) {
                    $count++;
                    @check[$rand]=0;
                    $seq = $1;
                    open (OUT1, ">>sample_seq_m.txt");
                    print (OUT1 ">sample_number$Z.cds\n");
                    close (OUT1);
                    
                    @str2 = split //, $seq;
                    
                    open (OUT1, ">>sample_seq_m.txt");
                    foreach $_ (@str2) {
                        if ($_ == 0) {print (OUT1 "T");}
                        if ($_ == 1) {print (OUT1 "C");}
                        if ($_ == 2) {print (OUT1 "A");}
                        if ($_ == 3) {print (OUT1 "G");}
                    }
                    print (OUT1 "\n");
                    close (OUT1);
                    $Z++;
                }
            }
        }
    }
}
open (IN1, "result_s.txt");
@str1 = <IN1>;
$count=1;
$Z=1;
for ($n1=0;$n1<=1000;$n1++) {
    @check[$n1]=1;
}
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /whole\ssequence/) {
        while ($count<=10) {
            $rand = int (rand 1000);
            if (@check[$rand]==1) {
                if (@str1[$n2+$rand+1] =~ /\=\s\"(\d\d\d\d\d+)/) {
                    $count++;
                    @check[$rand]=0;
                    $seq = $1;
                    open (OUT1, ">>sample_seq_s.txt");
                    print (OUT1 ">sample_number$Z.cds\n");
                    close (OUT1);
                    
                    @str2 = split //, $seq;
                    
                    open (OUT1, ">>sample_seq_s.txt");
                    foreach $_ (@str2) {
                        if ($_ == 0) {print (OUT1 "T");}
                        if ($_ == 1) {print (OUT1 "C");}
                        if ($_ == 2) {print (OUT1 "A");}
                        if ($_ == 3) {print (OUT1 "G");}
                    }
                    print (OUT1 "\n");
                    close (OUT1);
                    $Z++;
                }
            }
        }
    }
}

open (IN1, "result_t.txt");
@str1 = <IN1>;
$count=1;
$Z=1;
for ($n1=0;$n1<=1000;$n1++) {
    @check[$n1]=1;
}
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /whole\ssequence/) {
        while ($count<=10) {
            $rand = int (rand 1000);
            if (@check[$rand]==1) {
                if (@str1[$n2+$rand+1] =~ /\=\s\"(\d\d\d\d\d+)/) {
                    $count++;
                    @check[$rand]=0;
                    $seq = $1;
                    open (OUT1, ">>sample_seq_t.txt");
                    print (OUT1 ">sample_number$Z.cds\n");
                    close (OUT1);
                    
                    @str2 = split //, $seq;
                    
                    open (OUT1, ">>sample_seq_t.txt");
                    foreach $_ (@str2) {
                        if ($_ == 0) {print (OUT1 "T");}
                        if ($_ == 1) {print (OUT1 "C");}
                        if ($_ == 2) {print (OUT1 "A");}
                        if ($_ == 3) {print (OUT1 "G");}
                    }
                    print (OUT1 "\n");
                    close (OUT1);
                    $Z++;
                }
            }
        }
    }
}
open (IN1, "result_y.txt");
@str1 = <IN1>;
$count=1;
$Z=1;
for ($n1=0;$n1<=1000;$n1++) {
    @check[$n1]=1;
}
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /whole\ssequence/) {
        while ($count<=10) {
            $rand = int (rand 1000);
            if (@check[$rand]==1) {
                if (@str1[$n2+$rand+1] =~ /\=\s\"(\d\d\d\d\d+)/) {
                    $count++;
                    @check[$rand]=0;
                    $seq = $1;
                    open (OUT1, ">>sample_seq_y.txt");
                    print (OUT1 ">sample_number$Z.cds\n");
                    close (OUT1);
                    
                    @str2 = split //, $seq;
                    
                    open (OUT1, ">>sample_seq_y.txt");
                    foreach $_ (@str2) {
                        if ($_ == 0) {print (OUT1 "T");}
                        if ($_ == 1) {print (OUT1 "C");}
                        if ($_ == 2) {print (OUT1 "A");}
                        if ($_ == 3) {print (OUT1 "G");}
                    }
                    print (OUT1 "\n");
                    close (OUT1);
                    $Z++;
                }
            }
        }
    }
}
open (IN1, "result_e.txt");
@str1 = <IN1>;
$count=1;
$Z=1;
for ($n1=0;$n1<=1000;$n1++) {
    @check[$n1]=1;
}
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /whole\ssequence/) {
        while ($count<=10) {
            $rand = int (rand 1000);
            if (@check[$rand]==1) {
                if (@str1[$n2+$rand+1] =~ /\=\s\"(\d\d\d\d\d+)/) {
                    $count++;
                    @check[$rand]=0;
                    $seq = $1;
                    open (OUT1, ">>sample_seq_e.txt");
                    print (OUT1 ">sample_number$Z.cds\n");
                    close (OUT1);
                    
                    @str2 = split //, $seq;
                    
                    open (OUT1, ">>sample_seq_e.txt");
                    foreach $_ (@str2) {
                        if ($_ == 0) {print (OUT1 "T");}
                        if ($_ == 1) {print (OUT1 "C");}
                        if ($_ == 2) {print (OUT1 "A");}
                        if ($_ == 3) {print (OUT1 "G");}
                    }
                    print (OUT1 "\n");
                    close (OUT1);
                    $Z++;
                }
            }
        }
    }
}
open (IN1, "result_o.txt");
@str1 = <IN1>;
$count=1;
$Z=1;
for ($n1=0;$n1<=1000;$n1++) {
    @check[$n1]=1;
}
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /whole\ssequence/) {
        while ($count<=10) {
            $rand = int (rand 1000);
            if (@check[$rand]==1) {
                if (@str1[$n2+$rand+1] =~ /\=\s\"(\d\d\d\d\d+)/) {
                    $count++;
                    @check[$rand]=0;
                    $seq = $1;
                    open (OUT1, ">>sample_seq_o.txt");
                    print (OUT1 ">sample_number$Z.cds\n");
                    close (OUT1);
                    
                    @str2 = split //, $seq;
                    
                    open (OUT1, ">>sample_seq_o.txt");
                    foreach $_ (@str2) {
                        if ($_ == 0) {print (OUT1 "T");}
                        if ($_ == 1) {print (OUT1 "C");}
                        if ($_ == 2) {print (OUT1 "A");}
                        if ($_ == 3) {print (OUT1 "G");}
                    }
                    print (OUT1 "\n");
                    close (OUT1);
                    $Z++;
                }
            }
        }
    }
}




