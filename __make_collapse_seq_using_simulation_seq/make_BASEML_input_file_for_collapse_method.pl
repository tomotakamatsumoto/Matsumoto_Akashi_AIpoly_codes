open (IN1,"sample_seq.txt_collapse_m.MFA");
@str1 = <IN1>;
open (IN2,"sample_seq.txt_collapse_s.MFA");
@str2 = <IN2>;
open (IN3,"sample_seq.txt_collapse_t.MFA");
@str3 = <IN3>;
open (IN4,"sample_seq.txt_collapse_y.MFA");
@str4 = <IN4>;
open (IN5,"sample_seq.txt_collapse_e.MFA");
@str5 = <IN5>;
open (IN6,"sample_seq.txt_collapse_o.MFA");
@str6 = <IN6>;

open (IN7,"site_used_for_collapse_method.txt");
@str7 = <IN7>;

foreach $_ (@str1) {
    if ($_ =~ /collapse-1.cds/) {
        open (OUT1, ">>BASEML_input.MFA");
        printf (OUT1 ">mel_collapse-1.cds\n");
        close (OUT1);
    }
    if ($_ =~ /collapse-2.cds/) {
        open (OUT1, ">>BASEML_input.MFA");
        printf (OUT1 "\n>mel_collapse-2.cds\n");
        close (OUT1);
    }
    if ($_ =~ /^\w\w\w\w\w+/) {
        @STR1 = split //, $_;
        foreach $_ (@str7) {
            $num = $_;
            open (OUT1, ">>BASEML_input.MFA");
            printf (OUT1 "$STR1[$num]");
            close (OUT1);
        }
    }
}
foreach $_ (@str2) {
    if ($_ =~ /collapse-1.cds/) {
        open (OUT1, ">>BASEML_input.MFA");
        printf (OUT1 "\n>sim_collapse-1.cds\n");
        close (OUT1);
    }
    if ($_ =~ /collapse-2.cds/) {
        open (OUT1, ">>BASEML_input.MFA");
        printf (OUT1 "\n>sim_collapse-2.cds\n");
        close (OUT1);
    }
    if ($_ =~ /^\w\w\w\w\w+/) {
        @STR1 = split //, $_;
        foreach $_ (@str7) {
            $num = $_;
            open (OUT1, ">>BASEML_input.MFA");
            printf (OUT1 "$STR1[$num]");
            close (OUT1);
        }
    }
}
foreach $_ (@str3) {
    if ($_ =~ /collapse-1.cds/) {
        open (OUT1, ">>BASEML_input.MFA");
        printf (OUT1 "\n>tei_collapse-1.cds\n");
        close (OUT1);
    }
    if ($_ =~ /collapse-2.cds/) {
        open (OUT1, ">>BASEML_input.MFA");
        printf (OUT1 "\n>tei_collapse-2.cds\n");
        close (OUT1);
    }
    if ($_ =~ /^\w\w\w\w\w+/) {
        @STR1 = split //, $_;
        foreach $_ (@str7) {
            $num = $_;
            open (OUT1, ">>BASEML_input.MFA");
            printf (OUT1 "$STR1[$num]");
            close (OUT1);
        }
    }
}
foreach $_ (@str4) {
    if ($_ =~ /collapse-1.cds/) {
        open (OUT1, ">>BASEML_input.MFA");
        printf (OUT1 "\n>yak_collapse-1.cds\n");
        close (OUT1);
    }
    if ($_ =~ /collapse-2.cds/) {
        open (OUT1, ">>BASEML_input.MFA");
        printf (OUT1 "\n>yak_collapse-2.cds\n");
        close (OUT1);
    }
    if ($_ =~ /^\w\w\w\w\w+/) {
        @STR1 = split //, $_;
        foreach $_ (@str7) {
            $num = $_;
            open (OUT1, ">>BASEML_input.MFA");
            printf (OUT1 "$STR1[$num]");
            close (OUT1);
        }
    }
}
foreach $_ (@str5) {
    if ($_ =~ /collapse-1.cds/) {
        open (OUT1, ">>BASEML_input.MFA");
        printf (OUT1 "\n>ere_collapse-1.cds\n");
        close (OUT1);
    }
    if ($_ =~ /collapse-2.cds/) {
        open (OUT1, ">>BASEML_input.MFA");
        printf (OUT1 "\n>ere_collapse-2.cds\n");
        close (OUT1);
    }
    if ($_ =~ /^\w\w\w\w\w+/) {
        @STR1 = split //, $_;
        foreach $_ (@str7) {
            $num = $_;
            open (OUT1, ">>BASEML_input.MFA");
            printf (OUT1 "$STR1[$num]");
            close (OUT1);
        }
    }
}
foreach $_ (@str6) {
    if ($_ =~ /collapse-1.cds/) {
        open (OUT1, ">>BASEML_input.MFA");
        printf (OUT1 "\n>ore_collapse-1.cds\n");
        close (OUT1);
    }
    if ($_ =~ /collapse-2.cds/) {
        open (OUT1, ">>BASEML_input.MFA");
        printf (OUT1 "\n>ore_collapse-2.cds\n");
        close (OUT1);
    }
    if ($_ =~ /^\w\w\w\w\w+/) {
        @STR1 = split //, $_;
        foreach $_ (@str7) {
            $num = $_;
            open (OUT1, ">>BASEML_input.MFA");
            printf (OUT1 "$STR1[$num]");
            close (OUT1);
        }
    }
}
open (OUT1, ">>BASEML_input.MFA");
printf (OUT1 "\n");
close (OUT1);
