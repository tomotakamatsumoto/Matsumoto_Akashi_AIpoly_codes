# code to generate collapse sequence file of six species 


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
 
for((X=1;X<=100;X++))
 
 do
 
 perl translate_ATGC_10samples.pl;

 
 
 perl make_2mel_seq_for_collapse_method_nucleotide.pl;
 perl find_sites_with_less_than_three_states.pl;
 perl make_BASEML_input_file_for_collapse_method.pl;
 
 mkdir 100reps
 mkdir 100reps/"$X";
 mv sample_seq.txt_collapse_m.MFA 100reps/"$X"/sample_seq.txt_collapse_m.MFA;
 mv sample_seq.txt_collapse_s.MFA 100reps/"$X"/sample_seq.txt_collapse_s.MFA;
 mv sample_seq.txt_collapse_t.MFA 100reps/"$X"/sample_seq.txt_collapse_t.MFA;
 mv sample_seq.txt_collapse_y.MFA 100reps/"$X"/sample_seq.txt_collapse_y.MFA;
 mv sample_seq.txt_collapse_e.MFA 100reps/"$X"/sample_seq.txt_collapse_e.MFA;
 mv sample_seq.txt_collapse_o.MFA 100reps/"$X"/sample_seq.txt_collapse_o.MFA;
 
 mv site_with_more_than_three_states_m.txt 100reps/"$X"/site_with_more_than_three_states_m.txt;
 mv site_with_more_than_three_states_s.txt 100reps/"$X"/site_with_more_than_three_states_s.txt;
 mv site_with_more_than_three_states_t.txt 100reps/"$X"/site_with_more_than_three_states_t.txt;
 mv site_with_more_than_three_states_y.txt 100reps/"$X"/site_with_more_than_three_states_y.txt;
 mv site_with_more_than_three_states_e.txt 100reps/"$X"/site_with_more_than_three_states_e.txt;
 mv site_with_more_than_three_states_o.txt 100reps/"$X"/site_with_more_than_three_states_o.txt;
 
 mv sample_seq_m.txt 100reps/"$X"/sample_seq_m.txt;
 mv sample_seq_s.txt 100reps/"$X"/sample_seq_s.txt;
 mv sample_seq_t.txt 100reps/"$X"/sample_seq_t.txt;
 mv sample_seq_y.txt 100reps/"$X"/sample_seq_y.txt;
 mv sample_seq_e.txt 100reps/"$X"/sample_seq_e.txt;
 mv sample_seq_o.txt 100reps/"$X"/sample_seq_o.txt;
 
 mv site_used_for_collapse_method.txt 100reps/"$X"/site_used_for_collapse_method.txt;
 mv BASEML_input.MFA 100reps/"$X"/BASEML_input.MFA;
 
 
 mv site_not_used_for_collapse_method.txt 100reps/"$X"/site_not_used_for_collapse_method.txt;
 
 
done
 
rm result_m.txt;
rm result_s.txt;
rm result_t.txt;
rm result_y.txt;
rm result_e.txt;
rm result_o.txt;
rm result_ms.txt;
rm result_ty.txt;
rm result_eo.txt;
rm result_tyeo.txt;


