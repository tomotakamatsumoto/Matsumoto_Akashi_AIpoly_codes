open (IN1,"o_tmp.ctl");
@str1 = <IN1>;

foreach $_ (@str1) {
    if ($_ =~ /idum_init(\t+)=\s-15000/) {
        $R = int(rand 20000);
        open (OUT1, ">>o.ctl");
        print (OUT1 "idum_init					= -$R\n");
        close (OUT1);
    }
    else {
        open (OUT1, ">>o.ctl");
        print (OUT1 "$_");
        close (OUT1);
    }
}
open (IN1,"result.txt");
@str1 = <IN1>;

open (OUT1, ">>o.ctl");
print (OUT1 "\n");
close (OUT1);

foreach $_ (@str1) {
    if ($_ =~ /seq_\d/) {
        open (OUT1, ">>o.ctl");
        print (OUT1 "$_");
        close (OUT1);
    }
    if ($_ =~ /ancestor_/) {
        open (OUT1, ">>o.ctl");
        print (OUT1 "$_");
        close (OUT1);
    }
}



