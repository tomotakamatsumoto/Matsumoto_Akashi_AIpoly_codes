open (IN1,"site_with_more_than_three_states_m.txt");
@str1 = <IN1>;
open (IN2,"site_with_more_than_three_states_s.txt");
@str2 = <IN2>;
open (IN3,"site_with_more_than_three_states_t.txt");
@str3 = <IN3>;
open (IN4,"site_with_more_than_three_states_y.txt");
@str4 = <IN4>;
open (IN5,"site_with_more_than_three_states_e.txt");
@str5 = <IN5>;
open (IN6,"site_with_more_than_three_states_o.txt");
@str6 = <IN6>;



for ($n=0;$n<=99999;$n++) {
    $check=0;
    foreach $_ (@str1) {
        if ($_ == $n) {$check=1;}
    }
    foreach $_ (@str2) {
        if ($_ == $n) {$check=1;}
    }
    foreach $_ (@str3) {
        if ($_ == $n) {$check=1;}
    }
    foreach $_ (@str4) {
        if ($_ == $n) {$check=1;}
    }
    foreach $_ (@str5) {
        if ($_ == $n) {$check=1;}
    }
    foreach $_ (@str6) {
        if ($_ == $n) {$check=1;}
    }
    
    if ($check==0) {
        open (OUT1, ">>site_used_for_collapse_method.txt");
        printf (OUT1 "$n\n");
        close (OUT1);
    }
    elsif ($check==1) {
        open (OUT2, ">>site_not_used_for_collapse_method.txt");
        printf (OUT2 "$n\n");
        close (OUT2);
    }
}
