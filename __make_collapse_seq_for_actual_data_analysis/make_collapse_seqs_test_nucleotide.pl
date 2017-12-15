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
        $refpos1 = @{'SEQ'.${'st_group'.$n1}}[$x];
        $refpos2 = @{'SEQ'.${'st_group'.$n1}}[$x];

        
        $polys=0;
        for ($n2=${'st_group'.$n1}+1;$n2<=${'ed_group'.$n1};$n2++) {
            $pos = @{'SEQ'.$n2}[$x];
            
            if (($refpos1 ne $pos) && ($refpos2 ne $pos)) {
                $polys++;
                $refpos2 = $pos;
            }
        }
        if ($polys >= 2) {
            printf "error: more than two states at $xth polymorphic site\n";
        }
        
        if ($polys == 1) {
            
            $r = int (rand(10));

            if ($r<=4) {
                @{'nSEQ'.$collapse_index_A}[$x] = $refpos1;
                @{'nSEQ'.$collapse_index_B}[$x] = $refpos2;
            }
            elsif ($r>4) {
                @{'nSEQ'.$collapse_index_A}[$x] = $refpos2;
                @{'nSEQ'.$collapse_index_B}[$x] = $refpos1;
            }
        }
    
        if ($polys == 0) {
            
            @{'nSEQ'.$collapse_index_A}[$x] = $refpos1;
            @{'nSEQ'.$collapse_index_B}[$x] = $refpos1;
        }
    }
    $start = 2*$collapse_group_num+1;
    
    for ($n3=0;$n3<$single_seq_num;$n3++) {
        $ref= $start+$n3;
        @{'nSEQ'.$ref}[$x] = @{'SEQ'.${'single_seq'.$n3}}[$x];
    }

    
    $x++;
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


