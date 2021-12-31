import argparse
import numpy as np
import pandas as pd
import json
import networkx as nx
import math
from dateutil.rrule import rrule, DAILY
import datetime

def jhu_country_name_iso_match():
    """
    jhu github country name -> iso code
    :return:
    """
    country_name_to_iso = pd.read_csv(
        "https://raw.githubusercontent.com/AnthonyEbert/COVID-19_ISO-3166/master/JohnsHopkins-to-A3.csv")
    country_name_iso_match = {}
    for _, c in country_name_to_iso.iterrows():
        country_name_iso_match[c['Country/Region']] = c['alpha3']
    country_name_iso_match['Micronesia'] = 'FSM'
    with open('raw_data/jhu_country_name_iso_match.json', 'w') as f:
        json.dump(country_name_iso_match, f)


def jhu_raw_data_process(t0_date_jhu, data_pd):
    """
    :param t0_date_jhu: form 4/15/21
    :param data_pd: from jhu github
    :return:
    """
    with open('raw_data/jhu_country_name_iso_match.json', 'r') as f:
        country_name_iso_match = json.load(f)
    final_data = {}
    finished_c = []
    for _, data in data_pd.iterrows():
        d = data[t0_date_jhu]
        if not np.isnan(d):
            country_name = data['Country/Region']
            if country_name in country_name_iso_match:
                country_code = country_name_iso_match[country_name]
                if country_code != 'CHN':
                    if country_code not in final_data:
                        final_data[country_code] = d
                        if pd.isnull(data['Province/State']):
                            finished_c.append(country_code)
                    else:
                        if country_code not in finished_c:
                            final_data[country_code] += d
    return final_data


def jhu_process_china(t0_date_jhu, data_pd, data_dict):
    """
    t0_date_jhu = '4/15/21'
    data_pd = pd.read_csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data'
                          '/csse_covid_19_time_series/time_series_covid19_recovered_global.csv')
    data_dict = {}
    :param t0_date_jhu: time
    :param data_pd: full
    :param data_dict: dict
    :return:
    """
    data_greater_china = data_pd[data_pd['Country/Region'] == 'China']
    total_d = 0
    for _, data in data_greater_china.iterrows():
        d = data[t0_date_jhu]
        province = data['Province/State']
        if not pd.isnull(d) and not pd.isnull(province):
            if province == 'Hong Kong':
                data_dict['HKG'] = d
            elif province == 'Macau':
                data_dict['MAC'] = d
            else:
                total_d += d
    data_dict['CHN'] = total_d


def raw_jhu_data(t0_date_jhu):
    """
    jhu gihub -> cum confirmed cases, deaths, recovery number and cfr for each countries
    chn hkg mac additional
    :param t0_date_jhu:
    :return:
    """
    recovered_data = pd.read_csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data'
                                 '/csse_covid_19_time_series/time_series_covid19_recovered_global.csv')

    confirmed_data = pd.read_csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data'
                                 '/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv')

    death_data = pd.read_csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data'
                             '/csse_covid_19_time_series/time_series_covid19_deaths_global.csv')

    final_recovery_data = jhu_raw_data_process(t0_date_jhu, recovered_data)
    jhu_process_china(t0_date_jhu, recovered_data, final_recovery_data)
    final_confirmed_data = jhu_raw_data_process(t0_date_jhu, confirmed_data)
    jhu_process_china(t0_date_jhu, confirmed_data, final_confirmed_data)
    final_death_data = jhu_raw_data_process(t0_date_jhu, death_data)
    jhu_process_china(t0_date_jhu, death_data, final_death_data)

    cfr = {}
    for c in final_confirmed_data:
        if final_confirmed_data[c] != 0:
            cfr[c] = final_death_data[c] / final_confirmed_data[c]

    return final_recovery_data, final_confirmed_data, final_death_data, cfr


def raw_data_from_owid(t0_date_owid):
    """
    our world in data -> cum fully vaccinated people, population (year 2019), last_day, effective number
    :param t0_date_owid: example '2021-04-15'
    :return:
    """
    owid = pd.read_csv('https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv')
    owid = owid[owid['date'] == t0_date_owid]
    People_fully_vaccinated = {}
    Re = {}  # reproduction_rate
    Population = {}

    for _, data in owid.iterrows():
        country = data['iso_code']
        if len(country) == 3 and not np.isnan(data['population']):
            Population[country] = data['population']
            if not np.isnan(data['people_fully_vaccinated']):
                People_fully_vaccinated[country] = data['people_fully_vaccinated']
            elif not np.isnan(data['total_vaccinations']):
                People_fully_vaccinated[country] = math.floor(data['total_vaccinations'] / 2)
            else:
                People_fully_vaccinated[country] = 0
            Re[country] = data['reproduction_rate']
    return People_fully_vaccinated, Population, Re


def get_income_class_and_G_air():
    """
    gni->income high income + chn + rus
    G_air from OAG  (https://www.oag.com/)
    :return:
    dict for values: 1: HC 0: LMIC
    traffic graph sum mobility
    """
    with open('raw_data/GNI_info.json', 'r') as f:
        GNI_info = json.load(f)
    Income_class = {}

    for c in GNI_info:
        if GNI_info[c] >= 12536:
            Income_class[c] = 1
        else:
            Income_class[c] = 0
    Income_class['CHN'] = 1
    Income_class['RUS'] = 1
    G_air = nx.read_gexf('raw_data/G_air_2020.gexf')  # yearly data (sample data here)
    return Income_class, G_air



def get_common_countries(t0_date_jhu, t0_date_owid, t0_date_owid_last_day):
    """
    all cover
    :return:
    """

    # ---- raw data ------
    raw_income, raw_G_air = get_income_class_and_G_air()
    final_recovery_data, final_confirmed_data, final_death_data, cfr = raw_jhu_data(
        t0_date_jhu)

    People_fully_vaccinated, Population, Re = raw_data_from_owid(t0_date_owid)
    People_fully_vaccinated_last_day, _, _ = raw_data_from_owid(t0_date_owid_last_day)
    # --------------------

    c1 = set(list(raw_income.keys())) & set(list(People_fully_vaccinated.keys())) & set(list(Population.keys())) \
         & set(list(final_recovery_data.keys())) & set(list(final_confirmed_data.keys())) & set(
        list(final_death_data.keys())) & set(list(People_fully_vaccinated_last_day)) - set(
        list(['PER']))


    return raw_income, raw_G_air, final_recovery_data, final_confirmed_data, final_death_data, cfr, \
           People_fully_vaccinated, Population, Re, People_fully_vaccinated_last_day, list(c1)


def raw_dict_to_array(dict, common):
    """
    raw_dict -> array
    :param dict:
    :return:
    """
    array = np.ones(len(common)) * np.nan
    for c in dict:
        if c in common:
            array[common.index(c)] = dict[c]
    return array


def get_P(G, common):
    G.remove_edges_from(nx.selfloop_edges(G))  # remove self loop
    F = np.array(nx.adjacency_matrix(G, common, weight='weight').todense())
    F_temp = (np.tril(F, -1) + np.transpose(np.triu(F, 1))) / 2
    F = F_temp + F_temp.T
    Fn = np.sum(F, axis=1)
    P = np.true_divide(F, np.transpose(np.tile(Fn, (len(F), 1))))
    return P


def raw_graph_to_P(nxG, common):
    update_G = nx.DiGraph()
    for e in nxG.edges():
        if e[0] in common and e[1] in common:
            update_G.add_edge(e[0], e[1], weight=nxG[e[0]][e[1]]['weight'])
    return get_P(update_G, common)


def init_country_info(t0_date_jhu, t0_date_owid, t0_date_owid_last_day):
    raw_income, raw_G_air, final_recovery_data, final_confirmed_data, final_death_data, cfr, \
    People_fully_vaccinated, Population, Re, People_fully_vaccinated_last_day, common = \
        get_common_countries(t0_date_jhu, t0_date_owid, t0_date_owid_last_day)


    if 'CHN' in Re:
        Re['MAC'] = Re['CHN']
        Re['HKG'] = Re['CHN']
    Re_remain = raw_dict_to_array(Re, common)
    Re_remain_mean = np.nanmean(Re_remain)
    Re_remain[np.isnan(Re_remain) | (Re_remain == 0)] = Re_remain_mean

    country_info = {'N': raw_dict_to_array(Population, common),
                    'V': raw_dict_to_array(People_fully_vaccinated, common),
                    'R': raw_dict_to_array(final_recovery_data, common),
                    'D': raw_dict_to_array(final_death_data, common),
                    'IS': raw_dict_to_array(final_confirmed_data, common) -
                          raw_dict_to_array(final_recovery_data, common) - raw_dict_to_array(final_death_data, common),
                    'CFR': raw_dict_to_array(cfr, common),
                    'P': raw_graph_to_P(raw_G_air, common),
                    'Income': raw_dict_to_array(raw_income, common),
                    'common': common,
                    'Re': Re_remain,
                    'Vac_new': np.maximum((raw_dict_to_array(People_fully_vaccinated, common) -
                                           raw_dict_to_array(People_fully_vaccinated_last_day, common)) /
                                          raw_dict_to_array(Population, common), np.zeros(len(common)))}
    return country_info


def save_specific_date_info(t0_date_jhu, t0_date_owid, t0_date_owid_last_day):
    country_info = init_country_info(t0_date_jhu, t0_date_owid, t0_date_owid_last_day)
    country_info_new = {}
    for c in country_info:
        if c != 'P':
            country_info_new[c] = list(country_info[c])
        else:
            country_info_new[c] = country_info[c].tolist()
    path = 'specific_date_info/' + t0_date_owid + '.json'
    with open(path, 'w') as f:
        json.dump(country_info_new, f)


def read_specific_date_info(t0_date_owid):
    path = 'specific_date_info/' + t0_date_owid + '.json'
    with open(path, 'r') as f:
        country_info_temp = json.load(f)
    country_info = {c: np.array(country_info_temp[c]) for c in country_info_temp}
    return country_info


def get_max_vac_rate(model_start_date_owid):
    with open('covered_countries_iso.json', 'r') as f:
        all_c = json.load(f)
    need_time_series_rows_owid = [str(d)[:10] for d in
                                  list(rrule(freq=DAILY, dtstart=datetime.datetime.strptime('2020-01-01', '%Y-%m-%d'),
                                             until=datetime.datetime.strptime(model_start_date_owid, '%Y-%m-%d')))]
    owid_all = pd.read_csv(
        'https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv')
    owid_vaccination = owid_all.loc[
        owid_all['iso_code'].isin(all_c) & owid_all['date'].isin(need_time_series_rows_owid), ['iso_code', 'date',
                                                                                               'new_vaccinations_smoothed']]
    owid_vaccination_max = owid_vaccination.groupby('iso_code')['new_vaccinations_smoothed'].max().to_dict()
    with open('max_vaccinations.json', 'w') as f:
        json.dump(owid_vaccination_max, f)
# --------------------------------------------------------

def hyper():
    hyper_paras = {
        "t0_date_jhu": '6/15/21',  # model start date in jhu form
        "t0_date_owid": '2021-06-15',  # model start date in owid form
        "t0_date_owid_last_day": '2021-06-14',  # the date before the model start date in owid form [for calculate
        # vac rate]

        # source https://www.cdc.gov/coronavirus/2019-ncov/vaccines/fully-vaccinated.html
        "dv": 1,  # 14,  # the time needed for the body to build full immunity against the virus

        # https://science.sciencemag.org/content/370/6518/811.full
        "alpha": 1 / 5,  # 1/alpha is the infectious period

        # https://www.nature.com/articles/s41579-020-00459-7
        # "dd": 16,  # interval from onset to death

        # https://www.nature.com/articles/s41579-020-00459-7
        "sigma": 1 / 5,  # 1/\sigma is the incubation period

        # https://www.sciencedirect.com/science/article/pii/S2666524721000690?via%3Dihub
        "eta_base": 0.828,  # the base vaccine efficacy against infections

        # https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)00947-8/fulltext
        "epsilon_base": 0.967,  # the base vaccine efficacy against deaths

        # https://science.sciencemag.org/content/372/6540/363.full unclear 0.5 year / 1 year
        "varepsilon": 1 / 365,  # 1/\varepsilon is the duration of vaccinal immunity [days]

        "R0_basic": 5  # basic R0

    }
    return hyper_paras
