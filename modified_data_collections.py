import json
import hyper_paras_setting
import pandas as pd
import datetime
import numpy as np
import heapq
# modify data of LMICs
# input ------------------------------------
with open('input.json', 'r') as f:
    input_info = json.load(f)
hyper_paras = hyper_paras_setting.hyper()
# ------------------------------------------

# date_owid_to_matlab
def date_owid_to_matlab(owid_date):
    # transfer owid date format （'2021-06-15'）to the data format in matlab files ('15-Jun-2021')
    date = datetime.datetime.strptime(owid_date, '%Y-%m-%d')
    return date.strftime('%d-%b-%Y')


def modify_data(attr_before, attr_after, country_info, modified_info):
    # change R D I IFR in country_info
    # attr_before = 'R' # key in country_info
    # attr_after = 'corrected_R' # key in matlab csv
    for i, c in enumerate(country_info['common']):
        if c in modified_info[attr_after]:
            country_info[attr_before][i] = modified_info[attr_after][c]
    return country_info

def get_age_structure_info():
    # source path https://population.un.org/wpp/Download/Standard/Population/
    total_age_structure = [str(i * 5) + '-' + str((i + 1) * 5 - 1) for i in range(20)]
    total_age_structure.extend(['100+'])

    with open('raw_data/country_digital_iso_map.json', 'r') as f:
        country_digital_iso_map = json.load(f)

    age_structure_result = {}  #dict
    df = pd.read_excel('raw_data/age_strcuture_2020.xlsx')
    for row in df.iterrows():
        digital_code = str(row[1]['Country code'])
        if digital_code in country_digital_iso_map:
            country_code = country_digital_iso_map[digital_code]
            age_info = []
            for age_structure in total_age_structure:
                age_info.append(row[1][age_structure] * 1000)
            age_structure_result[country_code] = age_info
    return age_structure_result


def get_aggregate_age_info():
    age_structure_result = get_age_structure_info()
    aggregate_age_info = {}
    for c in age_structure_result:
        total_info = age_structure_result[c]
        aggregate_info = np.array([total_info[0],
                                   sum(total_info[1:4]),
                                   sum(total_info[4:6]),
                                   sum(total_info[6:8]),
                                   sum(total_info[8:10]),
                                   sum(total_info[10:13]),
                                   sum(total_info[13:15]),
                                   sum(total_info[15:17]),
                                   sum(total_info[17:])]) / np.sum(total_info)
        aggregate_age_info[c] = aggregate_info# * age_specific_IFR)
    return aggregate_age_info

def modify_country_info(model_start_date):
    # model_start_date = '2021-06-15'
    # vac_rate_limit_setting:
    # two settings ->  1. based on original max daily vac rate for LMICs and HICs
    #                  2. based on settings in Science paper https://www.science.org/doi/pdf/10.1126/science.abd7343 [weekly vaccination rate of 1% ]

    country_info = hyper_paras_setting.read_specific_date_info(model_start_date)
    country_info['c_init'] = country_info['Re'] * country_info['N'] / hyper_paras['R0_basic'] / (country_info['N']-country_info['IS']-country_info['R']-country_info['D'])
    del country_info['Re']

    modified_info = pd.read_csv('matlab_files/'+date_owid_to_matlab(model_start_date)
                                     + '.csv', delimiter=' ').set_index('iso_code').to_dict()

    # ------------ change R D active IFR for LMICs ---------------
    modify_data('R', 'corrected_R', country_info, modified_info)
    modify_data('D', 'corrected_D', country_info, modified_info)
    modify_data('IS', 'corrected_active', country_info, modified_info)
    modify_data('CFR', 'corrected_IFR', country_info, modified_info)
    # -----------------------------------------------------------

    # ------------ change CFR name ----------------------------
    country_info['CFR_old'] = country_info.pop('CFR')
    # --------------------------------------------------------

    # ----------- add max vac rate for LMICs & HICs -----------------
    with open('max_vaccinations.json', 'r') as f:
        owid_vaccination_max = json.load(f)
    HICs_vac_rate = []
    LMICs_vac_rate = []
    for c in owid_vaccination_max:
        index = list(country_info['common']).index(c)
        if country_info['Income'][index] == 1:
            HICs_vac_rate.append(owid_vaccination_max[c] / 2 / country_info['N'][index])
        else:
            LMICs_vac_rate.append(owid_vaccination_max[c] / 2 / country_info['N'][index])
    HIC_max = max(HICs_vac_rate)
    LMICs_max = heapq.nlargest(5, LMICs_vac_rate)[1]

    daily_vac_max = np.ones(len(country_info['N'])) * LMICs_max
    daily_vac_max[country_info['Income'] == 1] = HIC_max
    country_info['Vac_rate_max'] = daily_vac_max * country_info['N']
    # ---------------------------------------------------------------

    # ------------ add new CFR for LMICs & HICs ---------------------
    cfr_age_adjust = pd.read_excel('raw_data/CFR_age.xlsx')
    cfr_age_H = cfr_age_adjust['H'].values
    cfr_age_L = cfr_age_adjust['L'].values
    raw_age_structure = get_aggregate_age_info()
    CFR_new = []
    for i in range(len(country_info['common'])):
        c = country_info['common'][i]
        if c not in raw_age_structure:
            CFR_new.append(np.nan)
        else:
            if country_info['Income'][i] == 1:
                CFR_new.append(np.sum(raw_age_structure[c] * cfr_age_H))
            else:
                CFR_new.append(np.sum(raw_age_structure[c] * cfr_age_L))
    CFR_new = np.array(CFR_new)
    cfr_mean_H = np.nanmean(CFR_new[country_info['Income'] == 1])
    cfr_mean_L = np.nanmean(CFR_new[country_info['Income'] == 0])
    CFR_new[(country_info['Income'] == 1) & np.isnan(CFR_new)] = cfr_mean_H
    CFR_new[(country_info['Income'] == 0) & np.isnan(CFR_new)] = cfr_mean_L
    country_info['CFR'] = CFR_new
    # ---------------------------------------------------------------

    country_info_new = {}
    for c in country_info:
        if c != 'P':
            country_info_new[c] = list(country_info[c])
        else:
            country_info_new[c] = country_info[c].tolist()

    path = 'specific_date_info_modified/' + model_start_date + '.json'
    with open(path, 'w') as f:
        json.dump(country_info_new, f)
    return country_info

def read_specific_date_info(t0_date_owid):
    path = 'specific_date_info_modified/' + t0_date_owid + '.json'
    with open(path, 'r') as f:
        country_info_temp = json.load(f)
    country_info = {c: np.array(country_info_temp[c]) for c in country_info_temp}
    return country_info

with open('input.json', 'r') as f:
    input_info = json.load(f)

modify_country_info('2021-06-15')