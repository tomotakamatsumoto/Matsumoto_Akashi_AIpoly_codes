 
 perl 16_make_ancestral_site_probability_mstyeo_collapse_only_0%_prob_at_ms_weighted_and_AWP_filter_double_poly.pl;

 perl 16_make_estimated_site_frequency_spectrum_only_0%_prob_at_ms.pl;

 
 rm collapse.MFA.mfa;
 rm fltrst;
 rm anc_site_prob_m_0%_ms.txt;
 rm site_pattern_with_higher_than_0%_ms_node.txt;
 rm site_pattern_with_lower_than_0%_ms_node.txt;
 
  
 perl 16_make_actual_site_frequency_spectrum_collapse_only_0%_prob_at_ms.pl;

 
 rm sample_seq_m.txt;
 rm site_not_used_for_higher_than_0%_ms_analysis.txt;

 rm site_used_for_collapse_method.txt;
 rm BASEML_input_melpoly.MFA;
 
 perl 16_compare_actual_estimated_SFS.pl;
 
 rm actual_frequency_spectrum_m_0%_ms.txt;

 rm estimated_frequency_spectrum_m_0%_ms.txt;


