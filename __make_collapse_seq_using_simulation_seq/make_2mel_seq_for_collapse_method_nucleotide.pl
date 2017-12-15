open (IN1,"sample_seq_m.txt");
@str1 = <IN1>;
$x=1;
$remove_polyc=0;
$remove_polys=0;
$total_polyc=0;
$codnum=0;

foreach $_ (@str1) {

    if (($_ !~/\>/) && ($_ =~/^A||^T||^G||^C/)) {
        $seq[$x] = $_;
        $x++;
    }
}

for ($n=1;$n<=10;$n++) {
    @{'SEQ'.$n} = split //, $seq[$n];
}

$x=0;

while ($x<$#SEQ1) {
    $polyc=1;
    
    $ref1 = @{'SEQ'.1}[$x];
    $ref2 = @{'SEQ'.1}[$x];

    
    for ($n=2;$n<=10;$n++) {
        $pos = @{'SEQ'.$n}[$x];
        
        if (($ref1 ne $pos) && ($ref2 ne $pos)) {
            $polyc++;
            $ref2 = $pos;
        }
    }
    if ($polyc >= 2) {
        $total_polyc++;
        
        $r = int (rand(10));
        #printf "$r\n";
        if ($r<=4) {
            @{'nSEQ'.1}[$x] = $ref1;
            @{'nSEQ'.2}[$x] = $ref2;

        }
        elsif ($r>4) {
            @{'nSEQ'.1}[$x] = $ref2;
            @{'nSEQ'.2}[$x] = $ref1;
        }
    }
    
    if ($polyc == 1) {
        
        @{'nSEQ'.1}[$x] = $ref1;
        @{'nSEQ'.2}[$x] = $ref2;
    }
    if ($polyc > 2) {
        open (OUT2, ">>site_with_more_than_three_states_m.txt");
        printf (OUT2 "$x\n");
        close (OUT2)
    }
    $x++;
}
printf "$total_polyc\n";

open (OUT1, ">>sample_seq.txt_collapse_m.MFA");
print (OUT1 "\>collapse-1.cds\n");
$x=0;
while ($x<=$#nSEQ1) {
    if (@nSEQ3[$x] !~ /\s/) {
        print (OUT1 "@nSEQ1[$x]");
    }
    $x++;
}
print (OUT1 "\n");
print (OUT1 "\>collapse-2.cds\n");
$x=0;
while ($x<=$#nSEQ2) {
    if (@nSEQ4[$x] !~ /\s/) {
        print (OUT1 "@nSEQ2[$x]");
    }
    $x++;
}
print (OUT1 "\n");
close (OUT1);



open (IN1,"sample_seq_s.txt");
@str1 = <IN1>;
$x=1;
$remove_polyc=0;
$remove_polys=0;
$total_polyc=0;
$codnum=0;

foreach $_ (@str1) {
    
    if (($_ !~/\>/) && ($_ =~/^A||^T||^G||^C/)) {
        $seq[$x] = $_;
        $x++;
    }
}

for ($n=1;$n<=10;$n++) {
    @{'SEQ'.$n} = split //, $seq[$n];
}

$x=0;

while ($x<$#SEQ1) {
    $polyc=1;
    
    $ref1 = @{'SEQ'.1}[$x];
    $ref2 = @{'SEQ'.1}[$x];
    
    
    for ($n=2;$n<=10;$n++) {
        $pos = @{'SEQ'.$n}[$x];
        
        if (($ref1 ne $pos) && ($ref2 ne $pos)) {
            $polyc++;
            $ref2 = $pos;
        }
    }
    if ($polyc >= 2) {
        $total_polyc++;
        
        $r = int (rand(10));
        #printf "$r\n";
        if ($r<=4) {
            @{'nSEQ'.1}[$x] = $ref1;
            @{'nSEQ'.2}[$x] = $ref2;
            
        }
        elsif ($r>4) {
            @{'nSEQ'.1}[$x] = $ref2;
            @{'nSEQ'.2}[$x] = $ref1;
        }
    }
    
    if ($polyc == 1) {
        
        @{'nSEQ'.1}[$x] = $ref1;
        @{'nSEQ'.2}[$x] = $ref2;
    }
    if ($polyc > 2) {
        open (OUT2, ">>site_with_more_than_three_states_s.txt");
        printf (OUT2 "$x\n");
        close (OUT2)
    }
    $x++;
}
printf "$total_polyc\n";

open (OUT1, ">>sample_seq.txt_collapse_s.MFA");
print (OUT1 "\>collapse-1.cds\n");
$x=0;
while ($x<=$#nSEQ1) {
    if (@nSEQ3[$x] !~ /\s/) {
        print (OUT1 "@nSEQ1[$x]");
    }
    $x++;
}
print (OUT1 "\n");
print (OUT1 "\>collapse-2.cds\n");
$x=0;
while ($x<=$#nSEQ2) {
    if (@nSEQ4[$x] !~ /\s/) {
        print (OUT1 "@nSEQ2[$x]");
    }
    $x++;
}
print (OUT1 "\n");
close (OUT1);




open (IN1,"sample_seq_t.txt");
@str1 = <IN1>;
$x=1;
$remove_polyc=0;
$remove_polys=0;
$total_polyc=0;
$codnum=0;

foreach $_ (@str1) {
    
    if (($_ !~/\>/) && ($_ =~/^A||^T||^G||^C/)) {
        $seq[$x] = $_;
        $x++;
    }
}

for ($n=1;$n<=10;$n++) {
    @{'SEQ'.$n} = split //, $seq[$n];
}

$x=0;

while ($x<$#SEQ1) {
    $polyc=1;
    
    $ref1 = @{'SEQ'.1}[$x];
    $ref2 = @{'SEQ'.1}[$x];
    
    
    for ($n=2;$n<=10;$n++) {
        $pos = @{'SEQ'.$n}[$x];
        
        if (($ref1 ne $pos) && ($ref2 ne $pos)) {
            $polyc++;
            $ref2 = $pos;
        }
    }
    if ($polyc >= 2) {
        $total_polyc++;
        
        $r = int (rand(10));
        #printf "$r\n";
        if ($r<=4) {
            @{'nSEQ'.1}[$x] = $ref1;
            @{'nSEQ'.2}[$x] = $ref2;
            
        }
        elsif ($r>4) {
            @{'nSEQ'.1}[$x] = $ref2;
            @{'nSEQ'.2}[$x] = $ref1;
        }
    }
    
    if ($polyc == 1) {
        
        @{'nSEQ'.1}[$x] = $ref1;
        @{'nSEQ'.2}[$x] = $ref2;
    }
    if ($polyc > 2) {
        open (OUT2, ">>site_with_more_than_three_states_t.txt");
        printf (OUT2 "$x\n");
        close (OUT2)
    }
    $x++;
}
printf "$total_polyc\n";

open (OUT1, ">>sample_seq.txt_collapse_t.MFA");
print (OUT1 "\>collapse-1.cds\n");
$x=0;
while ($x<=$#nSEQ1) {
    if (@nSEQ3[$x] !~ /\s/) {
        print (OUT1 "@nSEQ1[$x]");
    }
    $x++;
}
print (OUT1 "\n");
print (OUT1 "\>collapse-2.cds\n");
$x=0;
while ($x<=$#nSEQ2) {
    if (@nSEQ4[$x] !~ /\s/) {
        print (OUT1 "@nSEQ2[$x]");
    }
    $x++;
}
print (OUT1 "\n");
close (OUT1);




open (IN1,"sample_seq_y.txt");
@str1 = <IN1>;
$x=1;
$remove_polyc=0;
$remove_polys=0;
$total_polyc=0;
$codnum=0;

foreach $_ (@str1) {
    
    if (($_ !~/\>/) && ($_ =~/^A||^T||^G||^C/)) {
        $seq[$x] = $_;
        $x++;
    }
}

for ($n=1;$n<=10;$n++) {
    @{'SEQ'.$n} = split //, $seq[$n];
}

$x=0;

while ($x<$#SEQ1) {
    $polyc=1;
    
    $ref1 = @{'SEQ'.1}[$x];
    $ref2 = @{'SEQ'.1}[$x];
    
    
    for ($n=2;$n<=10;$n++) {
        $pos = @{'SEQ'.$n}[$x];
        
        if (($ref1 ne $pos) && ($ref2 ne $pos)) {
            $polyc++;
            $ref2 = $pos;
        }
    }
    if ($polyc >= 2) {
        $total_polyc++;
        
        $r = int (rand(10));
        #printf "$r\n";
        if ($r<=4) {
            @{'nSEQ'.1}[$x] = $ref1;
            @{'nSEQ'.2}[$x] = $ref2;
            
        }
        elsif ($r>4) {
            @{'nSEQ'.1}[$x] = $ref2;
            @{'nSEQ'.2}[$x] = $ref1;
        }
    }
    
    if ($polyc == 1) {
        
        @{'nSEQ'.1}[$x] = $ref1;
        @{'nSEQ'.2}[$x] = $ref2;
    }
    if ($polyc > 2) {
        open (OUT2, ">>site_with_more_than_three_states_y.txt");
        printf (OUT2 "$x\n");
        close (OUT2)
    }
    $x++;
}
printf "$total_polyc\n";

open (OUT1, ">>sample_seq.txt_collapse_y.MFA");
print (OUT1 "\>collapse-1.cds\n");
$x=0;
while ($x<=$#nSEQ1) {
    if (@nSEQ3[$x] !~ /\s/) {
        print (OUT1 "@nSEQ1[$x]");
    }
    $x++;
}
print (OUT1 "\n");
print (OUT1 "\>collapse-2.cds\n");
$x=0;
while ($x<=$#nSEQ2) {
    if (@nSEQ4[$x] !~ /\s/) {
        print (OUT1 "@nSEQ2[$x]");
    }
    $x++;
}
print (OUT1 "\n");
close (OUT1);




open (IN1,"sample_seq_e.txt");
@str1 = <IN1>;
$x=1;
$remove_polyc=0;
$remove_polys=0;
$total_polyc=0;
$codnum=0;

foreach $_ (@str1) {
    
    if (($_ !~/\>/) && ($_ =~/^A||^T||^G||^C/)) {
        $seq[$x] = $_;
        $x++;
    }
}

for ($n=1;$n<=10;$n++) {
    @{'SEQ'.$n} = split //, $seq[$n];
}

$x=0;

while ($x<$#SEQ1) {
    $polyc=1;
    
    $ref1 = @{'SEQ'.1}[$x];
    $ref2 = @{'SEQ'.1}[$x];
    
    
    for ($n=2;$n<=10;$n++) {
        $pos = @{'SEQ'.$n}[$x];
        
        if (($ref1 ne $pos) && ($ref2 ne $pos)) {
            $polyc++;
            $ref2 = $pos;
        }
    }
    if ($polyc >= 2) {
        $total_polyc++;
        
        $r = int (rand(10));
        #printf "$r\n";
        if ($r<=4) {
            @{'nSEQ'.1}[$x] = $ref1;
            @{'nSEQ'.2}[$x] = $ref2;
            
        }
        elsif ($r>4) {
            @{'nSEQ'.1}[$x] = $ref2;
            @{'nSEQ'.2}[$x] = $ref1;
        }
    }
    
    if ($polyc == 1) {
        
        @{'nSEQ'.1}[$x] = $ref1;
        @{'nSEQ'.2}[$x] = $ref2;
    }
    if ($polyc > 2) {
        open (OUT2, ">>site_with_more_than_three_states_e.txt");
        printf (OUT2 "$x\n");
        close (OUT2)
    }
    $x++;
}
printf "$total_polyc\n";

open (OUT1, ">>sample_seq.txt_collapse_e.MFA");
print (OUT1 "\>collapse-1.cds\n");
$x=0;
while ($x<=$#nSEQ1) {
    if (@nSEQ3[$x] !~ /\s/) {
        print (OUT1 "@nSEQ1[$x]");
    }
    $x++;
}
print (OUT1 "\n");
print (OUT1 "\>collapse-2.cds\n");
$x=0;
while ($x<=$#nSEQ2) {
    if (@nSEQ4[$x] !~ /\s/) {
        print (OUT1 "@nSEQ2[$x]");
    }
    $x++;
}
print (OUT1 "\n");
close (OUT1);




open (IN1,"sample_seq_o.txt");
@str1 = <IN1>;
$x=1;
$remove_polyc=0;
$remove_polys=0;
$total_polyc=0;
$codnum=0;

foreach $_ (@str1) {
    
    if (($_ !~/\>/) && ($_ =~/^A||^T||^G||^C/)) {
        $seq[$x] = $_;
        $x++;
    }
}

for ($n=1;$n<=10;$n++) {
    @{'SEQ'.$n} = split //, $seq[$n];
}

$x=0;

while ($x<$#SEQ1) {
    $polyc=1;
    
    $ref1 = @{'SEQ'.1}[$x];
    $ref2 = @{'SEQ'.1}[$x];
    
    
    for ($n=2;$n<=10;$n++) {
        $pos = @{'SEQ'.$n}[$x];
        
        if (($ref1 ne $pos) && ($ref2 ne $pos)) {
            $polyc++;
            $ref2 = $pos;
        }
    }
    if ($polyc >= 2) {
        $total_polyc++;
        
        $r = int (rand(10));
        #printf "$r\n";
        if ($r<=4) {
            @{'nSEQ'.1}[$x] = $ref1;
            @{'nSEQ'.2}[$x] = $ref2;
            
        }
        elsif ($r>4) {
            @{'nSEQ'.1}[$x] = $ref2;
            @{'nSEQ'.2}[$x] = $ref1;
        }
    }
    
    if ($polyc == 1) {
        
        @{'nSEQ'.1}[$x] = $ref1;
        @{'nSEQ'.2}[$x] = $ref2;
    }
    if ($polyc > 2) {
        open (OUT2, ">>site_with_more_than_three_states_o.txt");
        printf (OUT2 "$x\n");
        close (OUT2)
    }
    $x++;
}
printf "$total_polyc\n";

open (OUT1, ">>sample_seq.txt_collapse_o.MFA");
print (OUT1 "\>collapse-1.cds\n");
$x=0;
while ($x<=$#nSEQ1) {
    if (@nSEQ3[$x] !~ /\s/) {
        print (OUT1 "@nSEQ1[$x]");
    }
    $x++;
}
print (OUT1 "\n");
print (OUT1 "\>collapse-2.cds\n");
$x=0;
while ($x<=$#nSEQ2) {
    if (@nSEQ4[$x] !~ /\s/) {
        print (OUT1 "@nSEQ2[$x]");
    }
    $x++;
}
print (OUT1 "\n");
close (OUT1);
