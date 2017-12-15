rm -r 10.input/AIpoly;

for X in 1


do

mkdir 10.input/AIpoly

cd 00.program

cp ../10.input/AIpoly_tmp/burn-in.ctl burn-in_tmp.ctl;
perl make_input_burn-in.pl;
mv burn-in.ctl ../10.input/AIpoly/burn-in.ctl;

./AIpoly_burnin < ../10.input/0_intput_MCU_AIpoly/input_burn-in;

cp ../20.output/fsim_out_1_035_04.txt result.txt;
mv ../20.output ../output_burn-in_"$X";
cp ../10.input/AIpoly_tmp/tyeo.ctl tyeo_tmp.ctl;
cp ../10.input/AIpoly_tmp/ms.ctl ms_tmp.ctl;

perl make_input_tyeo.pl;
mv tyeo.ctl ../10.input/AIpoly/tyeo.ctl;
perl make_input_ms.pl;
mv ms.ctl ../10.input/AIpoly/ms.ctl;
rm result.txt;

./AIpoly_after_burnin < ../10.input/0_intput_MCU_AIpoly/input_tyeo;

cp ../20.output/fsim_out_1_035_04.txt result.txt;
mv ../20.output ../output_tyeo_"$X";
cp ../10.input/AIpoly_tmp/ty.ctl ty_tmp.ctl;
cp ../10.input/AIpoly_tmp/eo.ctl eo_tmp.ctl;

perl make_input_ty.pl;
mv ty.ctl ../10.input/AIpoly/ty.ctl;
perl make_input_eo.pl;
mv eo.ctl ../10.input/AIpoly/eo.ctl;
rm result.txt;

./AIpoly_after_burnin < ../10.input/0_intput_MCU_AIpoly/input_ms;

cp ../20.output/fsim_out_1_035_04.txt result.txt;
mv ../20.output ../output_ms_"$X";
cp ../10.input/AIpoly_tmp/m.ctl m_tmp.ctl;
cp ../10.input/AIpoly_tmp/s.ctl s_tmp.ctl;

perl make_input_m.pl;
mv m.ctl ../10.input/AIpoly/m.ctl;
perl make_input_s.pl;
mv s.ctl ../10.input/AIpoly/s.ctl;
rm result.txt;

./AIpoly_after_burnin < ../10.input/0_intput_MCU_AIpoly/input_ty;

cp ../20.output/fsim_out_1_035_04.txt result.txt;
mv ../20.output ../output_ty_"$X";
cp ../10.input/AIpoly_tmp/t.ctl t_tmp.ctl;
cp ../10.input/AIpoly_tmp/y.ctl y_tmp.ctl;

perl make_input_t.pl;
mv t.ctl ../10.input/AIpoly/t.ctl;
perl make_input_y.pl;
mv y.ctl ../10.input/AIpoly/y.ctl;
rm result.txt;

./AIpoly_after_burnin < ../10.input/0_intput_MCU_AIpoly/input_eo;

cp ../20.output/fsim_out_1_035_04.txt result.txt;
mv ../20.output ../output_eo_"$X";
cp ../10.input/AIpoly_tmp/e.ctl e_tmp.ctl;
cp ../10.input/AIpoly_tmp/o.ctl o_tmp.ctl;

perl make_input_e.pl;
mv e.ctl ../10.input/AIpoly/e.ctl;
perl make_input_o.pl;
mv o.ctl ../10.input/AIpoly/o.ctl;
rm result.txt;

./AIpoly_after_burnin < ../10.input/0_intput_MCU_AIpoly/input_m;
mv ../20.output ../output_m_"$X";

./AIpoly_after_burnin < ../10.input/0_intput_MCU_AIpoly/input_s;
mv ../20.output ../output_s_"$X";

./AIpoly_after_burnin < ../10.input/0_intput_MCU_AIpoly/input_t;
mv ../20.output ../output_t_"$X";

./AIpoly_after_burnin < ../10.input/0_intput_MCU_AIpoly/input_y;
mv ../20.output ../output_y_"$X";

./AIpoly_after_burnin < ../10.input/0_intput_MCU_AIpoly/input_e;
mv ../20.output ../output_e_"$X";

./AIpoly_after_burnin < ../10.input/0_intput_MCU_AIpoly/input_o;
mv ../20.output ../output_o_"$X";

mv ../10.input/AIpoly ../10.input/AIpoly_"$X";

cd ..

done
