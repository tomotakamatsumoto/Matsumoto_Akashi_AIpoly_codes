
open (IN2,"ancestor_s.txt");
@str2 = <IN2>;
open (IN3,"ancestor_t.txt");
@str3 = <IN3>;
open (IN4,"ancestor_y.txt");
@str4 = <IN4>;
open (IN5,"ancestor_e.txt");
@str5 = <IN5>;
open (IN6,"ancestor_o.txt");
@str6 = <IN6>;
open (IN8,"ancestor_m.txt");
@str8 = <IN8>;

open (IN7,"site_used_for_collapse_method.txt");
@str7 = <IN7>;

open (OUT1, ">>BASEML_input_melpoly.MFA");

foreach $_ (@str8) {
    if ($_ =~ /^\d\d\d\d\d\d\d+/) {
        printf (OUT1 ">mel_anc.cds\n");
        @STR1 = split //, $_;
        foreach $_ (@str7) {
            $num = $_;
            if ($STR1[$num] == 0) {printf (OUT1 "T")};
            if ($STR1[$num] == 1) {printf (OUT1 "C")};
            if ($STR1[$num] == 2) {printf (OUT1 "A")};
            if ($STR1[$num] == 3) {printf (OUT1 "G")};
        }
    }
}
foreach $_ (@str2) {
    if ($_ =~ /^\d\d\d\d\d\d\d+/) {
        printf (OUT1 "\n>sim_anc.cds\n");
        @STR1 = split //, $_;
        foreach $_ (@str7) {
            $num = $_;
            if ($STR1[$num] == 0) {printf (OUT1 "T")};
            if ($STR1[$num] == 1) {printf (OUT1 "C")};
            if ($STR1[$num] == 2) {printf (OUT1 "A")};
            if ($STR1[$num] == 3) {printf (OUT1 "G")};
        }
    }
}
foreach $_ (@str3) {
    if ($_ =~ /^\d\d\d\d\d\d\d+/) {
        printf (OUT1 "\n>tei_anc.cds\n");
        @STR1 = split //, $_;
        foreach $_ (@str7) {
            $num = $_;
            if ($STR1[$num] == 0) {printf (OUT1 "T")};
            if ($STR1[$num] == 1) {printf (OUT1 "C")};
            if ($STR1[$num] == 2) {printf (OUT1 "A")};
            if ($STR1[$num] == 3) {printf (OUT1 "G")};
        }
    }
}
foreach $_ (@str4) {
    if ($_ =~ /^\d\d\d\d\d\d\d+/) {
        printf (OUT1 "\n>yak_anc.cds\n");
        @STR1 = split //, $_;
        foreach $_ (@str7) {
            $num = $_;
            if ($STR1[$num] == 0) {printf (OUT1 "T")};
            if ($STR1[$num] == 1) {printf (OUT1 "C")};
            if ($STR1[$num] == 2) {printf (OUT1 "A")};
            if ($STR1[$num] == 3) {printf (OUT1 "G")};
        }
    }
}
foreach $_ (@str5) {
    if ($_ =~ /^\d\d\d\d\d\d\d+/) {
        printf (OUT1 "\n>ere_anc.cds\n");
        @STR1 = split //, $_;
        foreach $_ (@str7) {
            $num = $_;
            if ($STR1[$num] == 0) {printf (OUT1 "T")};
            if ($STR1[$num] == 1) {printf (OUT1 "C")};
            if ($STR1[$num] == 2) {printf (OUT1 "A")};
            if ($STR1[$num] == 3) {printf (OUT1 "G")};
        }
    }
}
foreach $_ (@str6) {
    if ($_ =~ /^\d\d\d\d\d\d\d+/) {
        printf (OUT1 "\n>ore_anc.cds\n");
        @STR1 = split //, $_;
        foreach $_ (@str7) {
            $num = $_;
            if ($STR1[$num] == 0) {printf (OUT1 "T")};
            if ($STR1[$num] == 1) {printf (OUT1 "C")};
            if ($STR1[$num] == 2) {printf (OUT1 "A")};
            if ($STR1[$num] == 3) {printf (OUT1 "G")};
        }
    }
}
printf (OUT1 "\n");
close (OUT1);
