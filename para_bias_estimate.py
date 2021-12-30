import pandas as pd
import numpy as np
import datetime
import json
from dateutil.relativedelta import relativedelta



import hyper_paras_setting


def define_R_cal_start_date(model_start_date, months, days):
    """
    model_start_date: model start date(dateframe)
    months + days: "months" months + "days" days before
    return: dateframe
    """
    return model_start_date - relativedelta(months=months) - relativedelta(days=days)


def trans_string_to_time_jhu(date_jhu):
    """
    jhu date formate -> time for comparison
    :param date_jhu: example: 5/8/21
    :return: dateframe
    """
    mon_s, day_s, year_s = date_jhu.split('/')
    year_s = '20' + year_s
    mon_s = mon_s if len(mon_s) == 2 else '0' + mon_s
    day_s = day_s if len(day_s) == 2 else '0' + day_s
    return datetime.datetime(int(year_s), int(mon_s), int(day_s))


def trans_time_to_string_jhu(dateframe):
    """
    dateframe -> jhu date formate
    """
    year_s, mon_s, day_s = str(dateframe)[:10].split('-')
    year_s = year_s[-2:]
    mon_s = mon_s if mon_s[0] != '0' else mon_s[-1]
    day_s = day_s if day_s[0] != '0' else day_s[-1]
    return mon_s + '/' + day_s + '/' + year_s


def owid_date_to_jhu(date_owid):
    # date_owid = '2021-05-08'
    year, mon, day = date_owid.split('-')
    mon_jhu = mon if mon[0] != '0' else mon[1]
    day_jhu = day if day[0] != '0' else day[1]
    date_jhu = mon_jhu + '/' + day_jhu + '/' + year[-2:]
    return date_jhu


def date_series_data_jhu(raw_jhu_df, need_time_series_columns, need_country_iso):
    with open('raw_data/jhu_country_name_iso_match.json', 'r') as f:
        country_name_iso_match = json.load(f)
    need_columns = ['Country/Region']
    need_columns.extend(need_time_series_columns)
    df_update = pd.DataFrame(raw_jhu_df, columns=need_columns).fillna(0)
    df_sum = df_update.groupby(by='Country/Region').sum().reset_index().replace(
        {'Country/Region': country_name_iso_match})
    df_sum = df_sum.drop(df_sum[df_sum['Country/Region'].map(len) != 3].index)
    df_sum = df_sum.loc[df_sum['Country/Region'].isin(need_country_iso), :]
    return df_sum.set_index('Country/Region')


def save_original_series_for_init_correction(model_start_date):
    """
    # model_start_date = '2021-06-15'
    """
    country_info = hyper_paras_setting.read_specific_date_info(model_start_date)
    ISO = country_info['common']
    ISO_HIC = country_info['common'][country_info['Income'] == 1]
    ISO_LMIC = country_info['common'][country_info['Income'] == 0]
    R = country_info['R']
    D = country_info['D']
    I = country_info['IS'] + R + D

    # get test data from owid
    owid = pd.read_csv('https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv')
    owid_new = owid.loc[
        owid['iso_code'].isin(ISO) & (owid['date'] == model_start_date), ['iso_code',
                                                                          'date',
                                                                          'total_tests',
                                                                          'total_tests_per_thousand',
                                                                          'population']]

    # average test rate for non nan
    average_test_rate_among_tests_HIC = owid_new.loc[owid['iso_code'].isin(ISO_HIC), ['total_tests_per_thousand']][
        'total_tests_per_thousand'].mean()
    average_test_rate_among_tests_LMIC = owid_new.loc[owid['iso_code'].isin(ISO_LMIC), ['total_tests_per_thousand']][
        'total_tests_per_thousand'].mean()

    raw_R = []
    raw_I = []
    raw_D = []
    estimated_tests = []
    if_estimate_tests = []

    for _, row in owid_new.iterrows():
        c_index = list(ISO).index(row['iso_code'])
        confirmed = I[c_index]
        death = D[c_index]

        # Correct recovery data for several countries  & countries confirmed cases < deaths + recoveries
        recovery = R[c_index]

        if row['iso_code'] in ['THA', 'CMR', 'SRB',
                                 'USA', 'BEL', 'SWE', 'NLD', 'IRL', 'GBR', 'KNA',
                                 'GRC', 'ESP', 'FRA', 'FIN', 'NOR', 'CHE']:
            with open('raw_data_for_ini_corrections/' + row['iso_code'] + '_active.json', 'r') as f:
                c_active = int(json.load(f)[model_start_date])
                recovery = confirmed - c_active - death

        recovery = min(recovery, confirmed - death)


        if pd.isna(row['total_tests']):
            if country_info['Income'][c_index] == 1:
                estimated_tests.append(average_test_rate_among_tests_HIC * row['population'] / 1000)
            else:
                estimated_tests.append(average_test_rate_among_tests_LMIC * row['population'] / 1000)
            if_estimate_tests.append(True)
        else:
            estimated_tests.append(row['total_tests'])
            if_estimate_tests.append(False)

        raw_I.append(confirmed)
        raw_D.append(death)
        raw_R.append(recovery)

    owid_new['total recoveries'] = raw_R
    owid_new['total cases'] = raw_I
    owid_new['total deaths'] = raw_D
    owid_new['estimated_tests'] = estimated_tests
    owid_new['if_estimate_tests'] = if_estimate_tests

    owid_new.to_csv(
        'raw_data_for_ini_corrections/raw_' + str(model_start_date) + '.csv', index=False)

save_original_series_for_init_correction('2021-06-15')




