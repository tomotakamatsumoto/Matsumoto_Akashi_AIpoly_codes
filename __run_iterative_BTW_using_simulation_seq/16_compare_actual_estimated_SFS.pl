open (IN1,"actual_frequency_spectrum_m_0%_ms.txt");
@str1 = <IN1>;
open (IN2,"actual_frequency_spectrum_m_90%_ms.txt");
@str2 = <IN2>;
open (IN3,"actual_frequency_spectrum_m_95%_ms.txt");
@str3 = <IN3>;
open (IN4,"estimated_frequency_spectrum_m_0%_ms.txt");
@str4 = <IN4>;
open (IN5,"estimated_frequency_spectrum_m_90%_ms.txt");
@str5 = <IN5>;
open (IN6,"estimated_frequency_spectrum_m_95%_ms.txt");
@str6 = <IN6>;

$n=0;
foreach $_ (@str1) {
    if ($_ =~ /(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\n/) {
        @act_pu[$n] = $4+$5+$10+$12;
        @act_up[$n] = $1+$3+$8+$9;
        @act_n[$n] = $2+$6+$7+$11;
        $n++;
    }
}
$n=0;
foreach $_ (@str4) {
    if ($_ =~ /(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\t(\d.*)\n/) {
        @est_pu[$n] = $4+$5+$10+$12;
        @est_up[$n] = $1+$3+$8+$9;
        @est_n[$n] = $2+$6+$7+$11;
        $n++;
    }
}
open (OUT1, ">>est_SFS_m_0%_pu.txt");
open (OUT2, ">>est_SFS_m_0%_up.txt");
open (OUT3, ">>est_SFS_m_0%_n.txt");
open (OUT4, ">>act_SFS_m_0%_pu.txt");
open (OUT5, ">>act_SFS_m_0%_up.txt");
open (OUT6, ">>act_SFS_m_0%_n.txt");

$total_est_pu = 0;
$total_est_up = 0;
$total_est_n = 0;
$total_act_pu = 0;
$total_act_up = 0;
$total_act_n = 0;

for ($n=0;$n<=8;$n++) {
    $X = @est_pu[$n]-@act_pu[$n];
    $Y = @est_up[$n]-@act_up[$n];
    $Z = @est_n[$n]-@act_n[$n];
    $total_est_pu = $total_est_pu + @est_pu[$n];
    $total_est_up = $total_est_up + @est_up[$n];
    $total_est_n = $total_est_n + @est_n[$n];
    $total_act_pu = $total_act_pu + @act_pu[$n];
    $total_act_up = $total_act_up + @act_up[$n];
    $total_act_n = $total_act_n + @act_n[$n];
    print (OUT1 "@est_pu[$n]\t");
    print (OUT2 "@est_up[$n]\t");
    print (OUT3 "@est_n[$n]\t");
    print (OUT4 "@act_pu[$n]\t");
    print (OUT5 "@act_up[$n]\t");
    print (OUT6 "@act_n[$n]\t");
}
print (OUT1 "\t$total_est_pu\n");
print (OUT2 "\t$total_est_up\n");
print (OUT3 "\t$total_est_n\n");
print (OUT4 "\t$total_act_pu\n");
print (OUT5 "\t$total_act_up\n");
print (OUT6 "\t$total_act_n\n");
close (OUT1);
close (OUT2);
close (OUT3);
close (OUT4);
close (OUT5);
close (OUT6);

