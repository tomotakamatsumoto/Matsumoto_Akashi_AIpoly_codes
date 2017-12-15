open (IN1,"data.ffn");
@str1 = <IN1>;
$count=1;
$remove_polyc=0;
$remove_polys=0;
$total_polyc_m=0;
$total_polyc_s=0;
$codnum=0;

$collapse_group_num = $ARGV[0];

$x=0;
for ($n=0;$n<$collapse_group_num;$n++) {
    ${'st_group'.$n} = $ARGV[2*$n+1];
    ${'ed_group'.$n} = $ARGV[2*$n+2];
    $x = $x+2;
}

$single_seq_num = $ARGV[$x+1];

for ($n=0;$n<$single_seq_num;$n++) {
    ${'single_seq'.$n} = $ARGV[$x+2+$n];
}


foreach $_ (@str1) {
    if ($_ =~/\>.*\n/) {
        $index[$count] = $_;
    }
    
    if (($_ !~/\>/) && ($_ =~/^A||^T||^G||^C||^N||^\./)) {
        $seq[$count] = $_;
        $count++;
    }
}

for ($n=1;$n<=$count-1;$n++) {
    @{'SEQ'.$n} = split //, $seq[$n];
}

$x=0;

while ($x<$#SEQ1) {

    for ($n1=0;$n1<$collapse_group_num;$n1++) {

        $collapse_index_A = 2*$n1+1;
        $collapse_index_B = 2*$n1+2;
        
        $polyc=1;
        $refpos11 = @{'SEQ'.${'st_group'.$n1}}[$x];
        $refpos12 = @{'SEQ'.${'st_group'.$n1}}[$x+1];
        $refpos13 = @{'SEQ'.${'st_group'.$n1}}[$x+2];
        
        $refpos21 = @{'SEQ'.${'st_group'.$n1}}[$x];
        $refpos22 = @{'SEQ'.${'st_group'.$n1}}[$x+1];
        $refpos23 = @{'SEQ'.${'st_group'.$n1}}[$x+2];
        
        $refcodon1 = "$refpos11$refpos12$refpos13";
        $refcodon2 = "$refpos21$refpos22$refpos23";
        
        
        for ($n2=${'st_group'.$n1}+1;$n2<=${'ed_group'.$n1};$n2++) {
            $pos1 = @{'SEQ'.$n2}[$x];
            $pos2 = @{'SEQ'.$n2}[$x+1];
            $pos3 = @{'SEQ'.$n2}[$x+2];
            $codon = "$pos1$pos2$pos3";
            
            printf "$n1\t$codon\n";
            
            if (($refcodon1 ne "$pos1$pos2$pos3") && ($refcodon2 ne "$pos1$pos2$pos3")) {
                $polyc++;
                $refcodon2 = $codon;
                $refpos21 = $pos1;
                $refpos22 = $pos2;
                $refpos23 = $pos3;
            }
        }
        if ($polyc >= 2) {
            $total_polyc_m++;
            $polys=0;
            if ($refpos11 ne $refpos21) {
                $polys++;
            }
            if ($refpos12 ne $refpos22) {
                $polys++;
            }
            if ($refpos13 ne $refpos23) {
                $polys++;
            }
            if ($polys >= 1) {
                
                $r = int (rand(10));
                #printf "$r\n";
                if ($r<=4) {
                    @{'nSEQ'.$collapse_index_A}[$x] = $refpos11;
                    @{'nSEQ'.$collapse_index_A}[$x+1] = $refpos12;
                    @{'nSEQ'.$collapse_index_A}[$x+2] = $refpos13;
                    
                    @{'nSEQ'.$collapse_index_B}[$x] = $refpos21;
                    @{'nSEQ'.$collapse_index_B}[$x+1] = $refpos22;
                    @{'nSEQ'.$collapse_index_B}[$x+2] = $refpos23;
                }
                elsif ($r>4) {
                    @{'nSEQ'.$collapse_index_A}[$x] = $refpos21;
                    @{'nSEQ'.$collapse_index_A}[$x+1] = $refpos22;
                    @{'nSEQ'.$collapse_index_A}[$x+2] = $refpos23;
                    
                    @{'nSEQ'.$collapse_index_B}[$x] = $refpos11;
                    @{'nSEQ'.$collapse_index_B}[$x+1] = $refpos12;
                    @{'nSEQ'.$collapse_index_B}[$x+2] = $refpos13;
                }
            }
            elsif ($polys == 0) {
                printf "error: no poly in polymorphic site\n";
            }
        }
        
        if ($polyc == 1) {
            
            @{'nSEQ'.$collapse_index_A}[$x] = $refpos11;
            @{'nSEQ'.$collapse_index_A}[$x+1] = $refpos12;
            @{'nSEQ'.$collapse_index_A}[$x+2] = $refpos13;
            
            @{'nSEQ'.$collapse_index_B}[$x] = $refpos11;
            @{'nSEQ'.$collapse_index_B}[$x+1] = $refpos12;
            @{'nSEQ'.$collapse_index_B}[$x+2] = $refpos13;
        }
    }
    $start = 2*$collapse_group_num+1;
    
    for ($n3=0;$n3<$single_seq_num;$n3++) {
        $ref= $start+$n3;
        @{'nSEQ'.$ref}[$x] = @{'SEQ'.${'single_seq'.$n3}}[$x];
        @{'nSEQ'.$ref}[$x+1] = @{'SEQ'.${'single_seq'.$n3}}[$x+1];
        @{'nSEQ'.$ref}[$x+2] = @{'SEQ'.${'single_seq'.$n3}}[$x+2];
    }

    
    $x=$x+3;
}
open (OUT1, ">>total_polyc.txt");
print (OUT1 "$total_polyc_m\t$total_polyc_s\n");
close (OUT1);

open (OUT1, ">>ndata.ffn");

$start = 2*$collapse_group_num+1;
for ($n3=0;$n3<$single_seq_num;$n3++) {
    
    print (OUT1 "\>single\_seq\_$n3\n");
    $ref = $start+$n3;
    $x=0;
    while ($x<=$#nSEQ1) {
        if (@{'nSEQ'.$ref} !~ /\s/) {
            print (OUT1 "@{'nSEQ'.$ref}[$x]");
        }
        $x++;
    }
    print (OUT1 "\n");
}


for ($n3=0;$n3<$collapse_group_num;$n3++) {
    
    print (OUT1 "\>collapse\_seq\_$n3\_A\n");
    $refA = 2*$n3+1;
    $x=0;
    while ($x<=$#nSEQ1) {
        if (@{'nSEQ'.$refA} !~ /\s/) {
            print (OUT1 "@{'nSEQ'.$refA}[$x]");
        }
        $x++;
    }
    print (OUT1 "\n");
    print (OUT1 "\>collapse\_seq\_$n3\_B\n");
    $refB = 2*$n3+2;
    $x=0;
    while ($x<=$#nSEQ1) {
        if (@{'nSEQ'.$refB} !~ /\s/) {
            print (OUT1 "@{'nSEQ'.$refB}[$x]");
        }
        $x++;
    }
    print (OUT1 "\n");
}

close (OUT1);


