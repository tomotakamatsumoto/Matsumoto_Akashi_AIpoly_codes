rm -r 10.input/AIpoly;

for X in 1


do

mkdir 10.input/AIpoly

cd 00.program



cp ../output_ms_1/fsim_out_1_035_04.txt result.txt;

cp ../10.input/AIpoly_tmp/m.ctl m_tmp.ctl;
cp ../10.input/AIpoly_tmp/s.ctl s_tmp.ctl;

perl make_input_m.pl;
mv m.ctl ../10.input/AIpoly/m.ctl;

rm result.txt;

./AIpoly_after_burnin < ../10.input/0_intput_MCU_AIpoly/input_m;
mv ../20.output ../output_m_"$X";



cd ..

done
