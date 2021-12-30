model_start_date = '2021-06-15';
raw_data_path = strcat('../raw_data_for_ini_corrections/raw_', model_start_date, '.csv');
T = readtable(raw_data_path);
[c_all,~,~] = unique(T.iso_code);

% % 获取prior
e_sample_results = get_trunc_beta(0.8, 0.4^2,0.65,1);
p_sample_results = get_trunc_beta(0.99995, 0.01^2,0.9998,1);
q_MS_sample_results = get_trunc_beta(0.9, 0.2^2,0.85,1);
q_MN_sample_results = get_trunc_beta(0.25, 0.1^2,0.2,0.4);
P_MS_unt_sample_results = get_trunc_beta(0.3, 0.2^2,0.1,0.5);

[corrected_R,corrected_D,corrected_active] = get_corrected_info_for_specific_date_modify(model_start_date,c_all,T,e_sample_results,p_sample_results,q_MS_sample_results,q_MN_sample_results,P_MS_unt_sample_results,'', 1);


                                                                            