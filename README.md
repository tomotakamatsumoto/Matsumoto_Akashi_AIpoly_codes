# Matsumoto_Akashi_AIpoly_codes

"__simulation_codes_and_control_files_for_each_scenario" folder contains the C code and control file to run the simulation in each evolutionary scenario.
run "run_AIpoly_simulation.sh" to start the simulation.

"__make_collapse_seq_using_simulation_seq" folder contains the codes to make collapse sequences for each of six species using the simulation output.
run "run_make_BASEML_input_collapse_mstyeo_100_reps.sh" to start the program.
This folder also contains codes to make ancestral sequences of six species using the simulation output.
run "run_make_BASEML_input_6_ancestors_100_reps.sh" to start the program.

"__run_iterative_BTW_using_simulation_seq" folder contains the codes to conduct BTW analysis.
Please manually set the SFS used for the weighting in "16_make_ancestral_site_probability_mstyeo_collapse_only_0%_prob_at_ms_weighted_and_AWP_filter_double_poly.pl".
run "run_compare_actual_estimted_SFS_100reps" to start the program.
For the detail including the required input files, please check "memo.docx".


"__make_collapse_seq_for_actual_data_analysis" and "__run_iterative_BTW_for_actual_data_analysis" contain codes for BTW analysis applicable to the actual data.
We are still developing them.
The codes will be applicable arbitrary number of species and polymorphic sequences, and iterative BTWest analysis will be done automatically.
