open (IN1, "result_m.txt");
@str1 = <IN1>;
$num=0;
$numT=0;
$numC=0;
$numA=0;
$numG=0;
open (OUT1, ">>ancestor_m.txt");
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /SiteData/) {
        $num=1;
    }
    
    if ($num == 1) {
        while (@str1[$n2] =~ /(\d):0.0000000000:/g) {
            print (OUT1 "$1");
            if ($1 == 0) {$numT++;}
            if ($1 == 1) {$numC++;}
            if ($1 == 2) {$numA++;}
            if ($1 == 3) {$numG++;}
        }
    }
}
print (OUT1 "\n$numT\n$numC\n$numA\n$numG\n");
close (OUT1);

open (IN1, "result_s.txt");
@str1 = <IN1>;
$num=0;
$numT=0;
$numC=0;
$numA=0;
$numG=0;
 open (OUT1, ">>ancestor_s.txt");
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /SiteData/) {
        $num=1;
    }
    
    if ($num == 1) {
        while (@str1[$n2] =~ /(\d):0.0000000000:/g) {
            print (OUT1 "$1");
            if ($1 == 0) {$numT++;}
            if ($1 == 1) {$numC++;}
            if ($1 == 2) {$numA++;}
            if ($1 == 3) {$numG++;}
        }
    }
}
print (OUT1 "\n$numT\n$numC\n$numA\n$numG\n");
close (OUT1);

open (IN1, "result_t.txt");
@str1 = <IN1>;
$num=0;
$numT=0;
$numC=0;
$numA=0;
$numG=0;
open (OUT1, ">>ancestor_t.txt");
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /SiteData/) {
        $num=1;
    }
    
    if ($num == 1) {
        while (@str1[$n2] =~ /(\d):0.0000000000:/g) {
            print (OUT1 "$1");
            if ($1 == 0) {$numT++;}
            if ($1 == 1) {$numC++;}
            if ($1 == 2) {$numA++;}
            if ($1 == 3) {$numG++;}
        }
    }
}
print (OUT1 "\n$numT\n$numC\n$numA\n$numG\n");
close (OUT1);

open (IN1, "result_y.txt");
@str1 = <IN1>;
$num=0;
$numT=0;
$numC=0;
$numA=0;
$numG=0;
open (OUT1, ">>ancestor_y.txt");
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /SiteData/) {
        $num=1;
    }
    
    if ($num == 1) {
        while (@str1[$n2] =~ /(\d):0.0000000000:/g) {
            print (OUT1 "$1");
            if ($1 == 0) {$numT++;}
            if ($1 == 1) {$numC++;}
            if ($1 == 2) {$numA++;}
            if ($1 == 3) {$numG++;}
        }
    }
}
print (OUT1 "\n$numT\n$numC\n$numA\n$numG\n");
close (OUT1);

open (IN1, "result_e.txt");
@str1 = <IN1>;
$num=0;
$numT=0;
$numC=0;
$numA=0;
$numG=0;
open (OUT1, ">>ancestor_e.txt");
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /SiteData/) {
        $num=1;
    }
    
    if ($num == 1) {
        while (@str1[$n2] =~ /(\d):0.0000000000:/g) {
            print (OUT1 "$1");
            if ($1 == 0) {$numT++;}
            if ($1 == 1) {$numC++;}
            if ($1 == 2) {$numA++;}
            if ($1 == 3) {$numG++;}
        }
    }
}
print (OUT1 "\n$numT\n$numC\n$numA\n$numG\n");
close (OUT1);

open (IN1, "result_o.txt");
@str1 = <IN1>;
$num=0;
$numT=0;
$numC=0;
$numA=0;
$numG=0;
open (OUT1, ">>ancestor_o.txt");
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /SiteData/) {
        $num=1;
    }
    
    if ($num == 1) {
        while (@str1[$n2] =~ /(\d):0.0000000000:/g) {
            print (OUT1 "$1");
            if ($1 == 0) {$numT++;}
            if ($1 == 1) {$numC++;}
            if ($1 == 2) {$numA++;}
            if ($1 == 3) {$numG++;}
        }
    }
}
print (OUT1 "\n$numT\n$numC\n$numA\n$numG\n");
close (OUT1);

open (IN1, "result_ms.txt");
@str1 = <IN1>;
$num=0;
$numT=0;
$numC=0;
$numA=0;
$numG=0;
open (OUT1, ">>ancestor_ms.txt");
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /SiteData/) {
        $num=1;
    }
    
    if ($num == 1) {
        while (@str1[$n2] =~ /(\d):0.0000000000:/g) {
            print (OUT1 "$1");
            if ($1 == 0) {$numT++;}
            if ($1 == 1) {$numC++;}
            if ($1 == 2) {$numA++;}
            if ($1 == 3) {$numG++;}
        }
    }
}
print (OUT1 "\n$numT\n$numC\n$numA\n$numG\n");
close (OUT1);

open (IN1, "result_ty.txt");
@str1 = <IN1>;
$num=0;
$numT=0;
$numC=0;
$numA=0;
$numG=0;
open (OUT1, ">>ancestor_ty.txt");
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /SiteData/) {
        $num=1;
    }
    
    if ($num == 1) {
        while (@str1[$n2] =~ /(\d):0.0000000000:/g) {
            print (OUT1 "$1");
            if ($1 == 0) {$numT++;}
            if ($1 == 1) {$numC++;}
            if ($1 == 2) {$numA++;}
            if ($1 == 3) {$numG++;}
        }
    }
}
print (OUT1 "\n$numT\n$numC\n$numA\n$numG\n");
close (OUT1);

open (IN1, "result_eo.txt");
@str1 = <IN1>;
$num=0;
$numT=0;
$numC=0;
$numA=0;
$numG=0;
open (OUT1, ">>ancestor_eo.txt");
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /SiteData/) {
        $num=1;
    }
    
    if ($num == 1) {
        while (@str1[$n2] =~ /(\d):0.0000000000:/g) {
            print (OUT1 "$1");
            if ($1 == 0) {$numT++;}
            if ($1 == 1) {$numC++;}
            if ($1 == 2) {$numA++;}
            if ($1 == 3) {$numG++;}
        }
    }
}
print (OUT1 "\n$numT\n$numC\n$numA\n$numG\n");
close (OUT1);

open (IN1, "result_tyeo.txt");
@str1 = <IN1>;
$num=0;
$numT=0;
$numC=0;
$numA=0;
$numG=0;
open (OUT1, ">>ancestor_tyeo.txt");
for ($n2=0;$n2<=$#str1;$n2++) {
    if (@str1[$n2] =~ /SiteData/) {
        $num=1;
    }
    
    if ($num == 1) {
        while (@str1[$n2] =~ /(\d):0.0000000000:/g) {
            print (OUT1 "$1");
            if ($1 == 0) {$numT++;}
            if ($1 == 1) {$numC++;}
            if ($1 == 2) {$numA++;}
            if ($1 == 3) {$numG++;}
        }
    }
}
print (OUT1 "\n$numT\n$numC\n$numA\n$numG\n");
close (OUT1);

