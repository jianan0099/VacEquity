##### reproduction of data correction
1 run raw_data_for_ini_corrections/active_data_compl to get complementary data for the following countries [no update in recovery in JHU data]:
    need_correct_name_list = ['thailand', 'cameroon','serbia',
                           'us', 'belgium', 'sweden', 'netherlands', 'ireland', 'uk', 'saint-kitts-and-nevis',
                          'greece', 'spain', 'france', 'finland', 'norway', 'switzerland']
    need_correct_ISO_list = ['THA', 'CMR','SRB',
                              'USA', 'BEL', 'SWE', 'NLD', 'IRL', 'GBR', 'KNA',
                             'GRC', 'ESP', 'FRA', 'FIN', 'NOR', 'CHE']
    output data file format: ISO_active.json
2 run para_bias_estimate.py
- save_original_series_for_init_correction('2021-06-15') to get raw_data_for_ini_corrections/raw_2021-06-15.json
3 run matlab_files/para_bias.m to get corrected I R D active
5 run matlab_files/supple_bias.m to get figures in the Supplementary Information
6 run modify_country_info in modified_data_collections.py to get modified country info in '2021-06-15.json'
