open (IN1,"data1.ffn");
@str1 = <IN1>;
open (IN2,"data2.ffn");
@str2 = <IN2>;

$gene = $ARGV[0];

for ($n1=0;$n1<=5;$n1++) {
    ${'SEQ'.$n1} = @str1[2*$n1+1];
    @{'SEQ'.$n1} = split //, ${'SEQ'.$n1};
}


for ($n1=14;$n1<=34;$n1++) {
    ${'mSEQ'.$n1} = @str2[2*$n1+1];
    @{'mSEQ'.$n1} = split //, ${'mSEQ'.$n1};
}

for ($n=0;$n<=$#SEQ1-1;$n++) {

    
    $site = $n;
    $ref1 = @mSEQ14[$site];
    $ref2 = 'nan';
    $nref1 = 1;
    for ($n2=15;$n2<=34;$n2++) {
        
        if (@{'mSEQ'.$n2}[$site] ne $ref1) {
            if ($ref2 eq 'nan') {
                $ref2 = @{'mSEQ'.$n2}[$site];
            }
            elsif (@{'mSEQ'.$n2}[$site] ne $ref2) {
                printf "error: more than two state at site $site: $ref1 $ref2 @{'mSEQ'.$n2}[$site]\n";
            }
        }
    }

    if ($ref2 ne 'nan') {
        $check_AA = 0;
        $check_BB = 0;
        if ($ref1 eq @{'SEQ'.4}[$site]) {
            if ($ref2 eq @{'SEQ'.5}[$site]) {
                $check_AA++;
            }
        }
        
        if ($ref1 eq @{'SEQ'.5}[$site]) {
            if ($ref2 eq @{'SEQ'.4}[$site]) {
                $check_BB++;
            }
        }

        if (($check_AA==0) && ($check_BB==0)) {
            open (OUT1, ">>error_$gene.txt");
            print (OUT1 "miss\t$site\t$ref1\t$ref2\t@{'SEQ'.4}[$site]\t@{'SEQ'.5}[$site]\n");
            close (OUT1);
        }
        #if (($check_AA>0) || ($check_BB>0)) {
        #    printf "correct\t$site\t$ref1\t$ref2\t@{'SEQ'.4}[$site]\t@{'SEQ'.5}[$site]\n";
        #}
    }
}

