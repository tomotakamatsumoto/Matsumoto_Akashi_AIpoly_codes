# code to make the ancestral sequence of six species from simulation record
 
cp ../simulation_code_and_output/stationary_MCU=0.5/output_m_1/fsim_out_1_035_04.txt result_m.txt;
cp ../simulation_code_and_output/stationary_MCU=0.5/output_s_1/fsim_out_1_035_04.txt result_s.txt;
cp ../simulation_code_and_output/stationary_MCU=0.5/output_t_1/fsim_out_1_035_04.txt result_t.txt;
cp ../simulation_code_and_output/stationary_MCU=0.5/output_y_1/fsim_out_1_035_04.txt result_y.txt;
cp ../simulation_code_and_output/stationary_MCU=0.5/output_e_1/fsim_out_1_035_04.txt result_e.txt;
cp ../simulation_code_and_output/stationary_MCU=0.5/output_o_1/fsim_out_1_035_04.txt result_o.txt;
cp ../simulation_code_and_output/stationary_MCU=0.5/output_ms_1/fsim_out_1_035_04.txt result_ms.txt;
cp ../simulation_code_and_output/stationary_MCU=0.5/output_ty_1/fsim_out_1_035_04.txt result_ty.txt;
cp ../simulation_code_and_output/stationary_MCU=0.5/output_eo_1/fsim_out_1_035_04.txt result_eo.txt;
cp ../simulation_code_and_output/stationary_MCU=0.5/output_tyeo_1/fsim_out_1_035_04.txt result_tyeo.txt;

perl make_ancestoral_sequence.pl;

mkdir 100reps_ancestor;

for((X=1;X<=100;X++))

do

 cp -f 100reps/"$X"/site_used_for_collapse_method.txt site_used_for_collapse_method.txt;
 
 perl make_BASEML_input_file_6_ancestors.pl;
 
 rm site_used_for_collapse_method.txt;

 cd 100reps_ancestor;
 mkdir ancestor_"$X";
 cd ..;

 mv BASEML_input_melpoly.MFA 100reps_ancestor/ancestor_"$X"/BASEML_input_melpoly.MFA;
 

done

 