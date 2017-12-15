open (IN1,"burn-in_tmp.ctl");
@str1 = <IN1>;

foreach $_ (@str1) {
    if ($_ =~ /idum_init(\t+)=\s-15000/) {
        $R = int(rand 20000);
        open (OUT1, ">>burn-in.ctl");
        print (OUT1 "idum_init					= -$R\n");
        close (OUT1);
    }
    else {
        open (OUT1, ">>burn-in.ctl");
        print (OUT1 "$_");
        close (OUT1);
    }
}




