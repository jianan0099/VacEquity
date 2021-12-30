# set model para
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import hyper_paras_setting
import paras
import json
import modified_data_collections as mdc
import time

# input ------------------------------------
with open('input.json', 'r') as f:
    input_info = json.load(f)
hyper_paras = hyper_paras_setting.hyper()
# ------------------------------------------

# ------ read modified country info --------------------------------------
country_info = mdc.read_specific_date_info(hyper_paras['t0_date_owid'])
# ------------------------------------------------------------------------
popu_sum_H = sum(country_info['N'][country_info['Income'] == 1])
popu_sum_L = sum(country_info['N'][country_info['Income'] == 0])
popu_sum_all = sum(country_info['N'])


def run_overall_results(M=5, mu=1e-4, theta=0.26, Lambda=5,
                        c_thre1=0.3, c_thre2=1.0, NPI_change=0, NPI_change_tau=365 * 1, NPI_change_strong=7,
                        tau=365 * 5, d=15, Tau=365 * 0.5):

    results = {}
    end_time = {}
    virus_prop = {}

    for VAS_each in [1, 6, 7, 4]:
        model_eq = paras.DiscreteSVEIRD(Tau, VAS_each, 0, VAS_each, M, d, mu, theta, K=[1e-5, 0.5, 0], tau=tau,
                                        c_thre1=c_thre1, c_thre2=c_thre2, NPI_change=NPI_change, NPI_change_tau=NPI_change_tau,
                                        NPI_change_strong=NPI_change_strong, Lambda=Lambda)
        model_eq.trans_SVEIRD()

        results[str(VAS_each) + 'eq' + 'H_frac'] = np.array(model_eq.trans_length(model_eq.active_H)) / popu_sum_H
        results[str(VAS_each) + 'eq' + 'L_frac'] = np.array(model_eq.trans_length(model_eq.active_L)) / popu_sum_L
        results[str(VAS_each) + 'eq' + 'H_D_frac'] = np.array(model_eq.trans_length1(model_eq.D_H)) / popu_sum_H
        results[str(VAS_each) + 'eq' + 'L_D_frac'] = np.array(model_eq.trans_length1(model_eq.D_L)) / popu_sum_L
        results[str(VAS_each) + 'eq' + 'H_to_L'] = np.array(model_eq.trans_length1(model_eq.H_to_H))
        results[str(VAS_each) + 'eq' + 'L_to_H'] = np.array(model_eq.trans_length1(model_eq.L_to_H))
        results[str(VAS_each) + 'eq' + 'HL_frac_H'] = np.array(model_eq.trans_length(model_eq.HL_frac_H))
        results[str(VAS_each) + 'eq' + 'HL_frac_L'] = np.array(model_eq.trans_length(model_eq.HL_frac_L))
        results[str(VAS_each) + 'eq' + 'cum'] = np.array(
            model_eq.trans_length1(model_eq.cum_H)) / popu_sum_all + np.array(
            model_eq.trans_length1(model_eq.cum_L)) / popu_sum_all
        results[str(VAS_each) + 'eq' + 'cum_D'] = np.array(
            model_eq.trans_length1(model_eq.D_H)) / popu_sum_all + np.array(
            model_eq.trans_length1(model_eq.D_L)) / popu_sum_all
        results[str(VAS_each) + 'eq' + 'Re_H'] = np.array(model_eq.trans_length(model_eq.Re_H))
        results[str(VAS_each) + 'eq' + 'Re_L'] = np.array(model_eq.trans_length(model_eq.Re_L))
        end_time[str(VAS_each) + 'eq'] = [model_eq.pandemic_end_time]
        for m in range(M):
            strain_frac = list(np.array(model_eq.cum1)[:, m])
            strain_frac.extend([0] * (tau - len(strain_frac)))
            results[str(VAS_each) + 'eq' + str(m + 1)] = np.array(strain_frac) / popu_sum_all
            results[str(VAS_each) + 'eq' + str(m + 1) + 'H'] = np.array(
                model_eq.trans_length(np.array(model_eq.strain_H)[:, m])) / popu_sum_H
            results[str(VAS_each) + 'eq' + str(m + 1) + 'L'] = np.array(
                model_eq.trans_length(np.array(model_eq.strain_L)[:, m])) / popu_sum_L

        for chi in [0.7, 0.8, 0.9]:
            # inequity
            model = paras.DiscreteSVEIRD(Tau, 100, chi, VAS_each, M, d, mu, theta, K=[1e-5, 0.5, 0], tau=tau,
                                         c_thre1=c_thre1, c_thre2=c_thre2, NPI_change=NPI_change, NPI_change_tau=NPI_change_tau,
                                         NPI_change_strong=NPI_change_strong, Lambda=Lambda)
            model.trans_SVEIRD()

            results[str(VAS_each) + 'ineq' + str(chi) + 'H_frac'] = np.array(
                model.trans_length(model.active_H)) / popu_sum_H
            results[str(VAS_each) + 'ineq' + str(chi) + 'L_frac'] = np.array(
                model.trans_length(model.active_L)) / popu_sum_L
            results[str(VAS_each) + 'ineq' + str(chi) + 'H_D_frac'] = np.array(
                model.trans_length1(model.D_H)) / popu_sum_H
            results[str(VAS_each) + 'ineq' + str(chi) + 'L_D_frac'] = np.array(
                model.trans_length1(model.D_L)) / popu_sum_L
            results[str(VAS_each) + 'ineq' + str(chi) + 'H_to_L'] = np.array(model.trans_length1(model.H_to_H))
            results[str(VAS_each) + 'ineq' + str(chi) + 'L_to_H'] = np.array(model.trans_length1(model.L_to_H))
            results[str(VAS_each) + 'ineq' + str(chi) + 'cum'] = np.array(
                model.trans_length1(model.cum_H)) / popu_sum_all + np.array(
                model.trans_length1(model.cum_L)) / popu_sum_all
            results[str(VAS_each) + 'ineq' + str(chi) + 'cum_D'] = np.array(
                model.trans_length1(model.D_H)) / popu_sum_all + np.array(model.trans_length1(model.D_L)) / popu_sum_all
            results[str(VAS_each) + 'ineq' + str(chi) + 'Re_H'] = np.array(model.trans_length(model.Re_H))
            results[str(VAS_each) + 'ineq' + str(chi) + 'Re_L'] = np.array(model.trans_length(model.Re_L))
            results[str(VAS_each) + 'ineq' + str(chi) + 'HL_frac_H'] = np.array(model.trans_length(model.HL_frac_H))
            results[str(VAS_each) + 'ineq' + str(chi) + 'HL_frac_L'] = np.array(model.trans_length(model.HL_frac_L))
            end_time[str(VAS_each) + 'ineq' + str(chi)] = [model.pandemic_end_time]
            for m in range(M):
                strain_frac = list(np.array(model.cum1)[:, m])
                strain_frac.extend([0] * (tau - len(strain_frac)))
                results[str(VAS_each) + 'ineq' + str(chi) + str(m + 1)] = np.array(strain_frac) / popu_sum_all
                results[str(VAS_each) + 'ineq' + str(chi) + str(m + 1) + 'H'] = np.array(
                    model.trans_length(np.array(model.strain_H)[:, m])) / popu_sum_H
                results[str(VAS_each) + 'ineq' + str(chi) + str(m + 1) + 'L'] = np.array(
                    model.trans_length(np.array(model.strain_L)[:, m])) / popu_sum_L
    df = pd.DataFrame(results)
    df.to_csv('matlab_files/results/overall_' + str(M) + '_' + str(mu) + '_' + str(theta) + '_' + str(Lambda) + '_' +
              str(d) + '_' +
              str(c_thre1) + '_' + str(c_thre2) + '_' + str(NPI_change) + '_' + str(NPI_change_tau) + '_' + str(NPI_change_strong) +
              '_' + input_info['vac_max_rate_choice'] + '_' + input_info['if_wanned_natural_immunity'] + '.csv',
              index=False)
    df = pd.DataFrame(end_time)
    df.to_csv('matlab_files/results/endtime_' + str(M) + '_' + str(mu) + '_' + str(theta) + '_' + str(Lambda) + '_' +
              str(d) + '_' +
              str(c_thre1) + '_' + str(c_thre2) + '_' + str(NPI_change) + '_' + str(NPI_change_tau) + '_' + str(NPI_change_strong) +
              '_' + input_info['vac_max_rate_choice'] + '_' + input_info['if_wanned_natural_immunity'] + '.csv',
              index=False)
    virus_prop['T'] = np.diagonal(model.T)
    virus_prop['F'] = np.mean(model.F, axis=0)
    virus_prop['eta'] = model.eta
    virus_prop['epsilon'] = model.epsilon
    df = pd.DataFrame(virus_prop)
    df.to_csv('matlab_files/results/virus_prop' + str(M) + '_' + str(theta) + '_'+ str(d) +'.csv')



def donation_compare(M=4, mu=5.6e-3, theta=0.26, Lambda=100,
                     Tau=365 * 0.5, d=15, tau=365 * 5, chi=0.8, VAS_each=1, VAS=100,
                     c_thre1=0.3, c_thre2=1.0, NPI_change=0, NPI_change_tau=365 * 1, NPI_change_strong=7):
    rich_num = sum(country_info['Income'] == 1)
    poor_num = sum(country_info['Income'] == 0)


    H_benefit_num = []
    H_benefit_sum = []
    L_benefit_num = []
    L_benefit_sum = []
    H_not_benefit_sum = []
    L_compare_to_eq_sum = []
    L_not_benefit_sum = []
    H_benefit_ave = []
    L_benefit_ave = []

    # ----- add deaths ------------------
    D_H_benefit_num = []
    D_H_benefit_sum = []
    D_L_benefit_num = []
    D_L_benefit_sum = []
    D_H_not_benefit_sum = []
    D_L_compare_to_eq_sum = []
    D_H_benefit_ave = []
    D_L_benefit_ave = []
    # -----------------------------------


    model_eq = paras.DiscreteSVEIRD(Tau, VAS_each, 0, VAS_each, M, d, mu, theta, K=[0, 0.5, 0], tau=tau,
                                    c_thre1=c_thre1, c_thre2=c_thre2, NPI_change=NPI_change, NPI_change_tau=NPI_change_tau,
                                    NPI_change_strong=NPI_change_strong, Lambda=Lambda)
    model_eq.trans_SVEIRD()
    cum_eq_L = model_eq.cum_all[country_info['Income'] == 0]
    cum_D_eq_L = model_eq.cum_D_all[country_info['Income'] == 0]
    # ------------------------------------------------------


    model_ineq = paras.DiscreteSVEIRD(Tau, VAS, chi, VAS_each, M, d, mu, theta, K=[0, 0.5, 0], tau=tau,
                                 c_thre1=c_thre1, c_thre2=c_thre2, NPI_change=NPI_change, NPI_change_tau=NPI_change_tau,
                                 NPI_change_strong=NPI_change_strong, Lambda=Lambda)
    model_ineq.trans_SVEIRD()
    cum_ineq_H = model_ineq.cum_all[country_info['Income'] == 1]
    cum_ineq_L = model_ineq.cum_all[country_info['Income'] == 0]
    cum_D_ineq_H = model_ineq.cum_D_all[country_info['Income'] == 1]
    cum_D_ineq_L = model_ineq.cum_D_all[country_info['Income'] == 0]
    # ------------------------------------------------------

    for K1 in np.arange(0, 1e-4 + 4e-8, 2e-6):
        H_benefit_num_temp = []
        H_benefit_sum_temp = []
        H_not_benefit_sum_temp = []
        L_benefit_num_temp = []
        L_benefit_sum_temp = []
        L_compare_to_eq_sum_temp = []
        L_not_benefit_sum_temp = []
        H_benefit_ave_temp = []
        L_benefit_ave_temp = []

        # ----- add death info -------------
        D_H_benefit_num_temp = []
        D_H_benefit_sum_temp = []
        D_H_not_benefit_sum_temp = []
        D_L_benefit_num_temp = []
        D_L_benefit_sum_temp = []
        D_L_compare_to_eq_sum_temp = []
        D_H_benefit_ave_temp = []
        D_L_benefit_ave_temp = []
        # --------------------------------

        for K2 in np.arange(0, 1.01, 0.02):
            aa = time.time()
            print(K1, K2)
            model_ineq_ad = paras.DiscreteSVEIRD(Tau, VAS, chi, VAS_each, M, d, mu, theta, K=[K1, K2, 1], tau=tau,
                                              c_thre1=c_thre1, c_thre2=c_thre2, NPI_change=NPI_change, NPI_change_tau=NPI_change_tau,
                                              NPI_change_strong=NPI_change_strong, Lambda=Lambda)

            model_ineq_ad.trans_SVEIRD()
            cum_ineq_ad_H = model_ineq_ad.cum_all[country_info['Income'] == 1]
            cum_ineq_ad_L = model_ineq_ad.cum_all[country_info['Income'] == 0]
            cum_D_ineq_ad_H = model_ineq_ad.cum_D_all[country_info['Income'] == 1]
            cum_D_ineq_ad_L = model_ineq_ad.cum_D_all[country_info['Income'] == 0]
            # --------------------------------------------------------------


            benefit_H = cum_ineq_H > cum_ineq_ad_H
            not_benefit_H = cum_ineq_H <= cum_ineq_ad_H
            benefit_L = cum_ineq_L > cum_ineq_ad_L
            not_benefit_L = cum_ineq_L <= cum_ineq_ad_L
            # -------------------------------------------------------------


            D_benefit_H = cum_D_ineq_H > cum_D_ineq_ad_H
            D_not_benefit_H = cum_D_ineq_H <= cum_D_ineq_ad_H
            D_benefit_L = cum_D_ineq_L > cum_D_ineq_ad_L
            # -------------------------------------------------------------


            H_benefit_num_temp.append(sum(benefit_H) / rich_num)  # cases
            D_H_benefit_num_temp.append(sum(D_benefit_H) / rich_num)  # deaths
            # -------------------------------------------------------------

            H_benefit_ave_temp.append(np.mean((cum_ineq_H-cum_ineq_ad_H)/country_info['N'][country_info['Income'] == 1]))
            L_benefit_ave_temp.append(
                np.mean((cum_ineq_L - cum_ineq_ad_L) / country_info['N'][country_info['Income'] == 0]))
            D_H_benefit_ave_temp.append(np.mean((cum_D_ineq_H-cum_D_ineq_ad_H)/country_info['N'][country_info['Income'] == 1]))
            D_L_benefit_ave_temp.append(
                np.mean((cum_D_ineq_L - cum_D_ineq_ad_L) / country_info['N'][country_info['Income'] == 0]))

            if sum(benefit_H) > 0:
                H_benefit_sum_temp.append(np.mean(
                    (cum_ineq_H[benefit_H] - cum_ineq_ad_H[benefit_H]) / country_info['N'][country_info['Income'] == 1][
                        benefit_H]))
            else:
                H_benefit_sum_temp.append(0)

            if sum(not_benefit_H) > 0:
                H_not_benefit_sum_temp.append(np.mean((cum_ineq_ad_H[not_benefit_H] - cum_ineq_H[not_benefit_H]) /
                                                      country_info['N'][country_info['Income'] == 1][
                                                          not_benefit_H]))
            else:
                H_not_benefit_sum_temp.append(0)

            # ------------ add deaths -----------------------------
            if sum(D_benefit_H) > 0:
                D_H_benefit_sum_temp.append(np.mean(
                    (cum_D_ineq_H[D_benefit_H] - cum_D_ineq_ad_H[D_benefit_H]) /
                    country_info['N'][country_info['Income'] == 1][D_benefit_H]))
            else:
                D_H_benefit_sum_temp.append(0)
            if sum(D_not_benefit_H) > 0:
                D_H_not_benefit_sum_temp.append(np.mean(
                    (cum_D_ineq_ad_H[D_not_benefit_H] - cum_D_ineq_H[D_not_benefit_H]) /
                    country_info['N'][country_info['Income'] == 1][
                        D_not_benefit_H]))
            else:
                D_H_not_benefit_sum_temp.append(0)
            # -----------------------------------------------------

            L_benefit_num_temp.append(sum(benefit_L) / poor_num)

            if sum(benefit_L) > 0:
                L_benefit_sum_temp.append(
                    np.mean((cum_ineq_L[benefit_L] - cum_ineq_ad_L[benefit_L]) /
                            country_info['N'][country_info['Income'] == 0][benefit_L]))
            else:
                L_benefit_sum_temp.append(0)


            L_compare_to_eq_sum_temp.append(
                np.mean(np.abs(cum_ineq_ad_L - cum_eq_L) / country_info['N'][country_info['Income'] == 0]))

            if sum(not_benefit_L) > 0:
                L_not_benefit_sum_temp.append(np.mean((cum_ineq_ad_L[not_benefit_L] - cum_ineq_L[not_benefit_L]) /
                                                      country_info['N'][country_info['Income'] == 0][not_benefit_L]))
            else:
                L_not_benefit_sum_temp.append(0)

            # ------------ add deaths -----------------------------
            D_L_benefit_num_temp.append(sum(D_benefit_L) / poor_num)
            if sum(D_benefit_L) > 0:
                D_L_benefit_sum_temp.append(
                    np.mean((cum_D_ineq_L[D_benefit_L] - cum_D_ineq_ad_L[D_benefit_L]) /
                            country_info['N'][country_info['Income'] == 0][D_benefit_L]))
            else:
                D_L_benefit_sum_temp.append(0)
            D_L_compare_to_eq_sum_temp.append(
                np.mean(np.abs((cum_D_ineq_ad_L - cum_D_eq_L)) / country_info['N'][country_info['Income'] == 0]))
            # -----------------------------------------------------
            print('time', time.time() - aa)
        H_benefit_num.append(H_benefit_num_temp)
        H_benefit_sum.append(H_benefit_sum_temp)
        L_benefit_num.append(L_benefit_num_temp)
        L_benefit_sum.append(L_benefit_sum_temp)
        H_not_benefit_sum.append(H_not_benefit_sum_temp)
        L_compare_to_eq_sum.append(L_compare_to_eq_sum_temp)
        L_not_benefit_sum.append(L_not_benefit_sum_temp)
        H_benefit_ave.append(H_benefit_ave_temp)
        L_benefit_ave.append(L_benefit_ave_temp)

        # ------------- add deaths --------------------
        D_H_benefit_num.append(D_H_benefit_num_temp)
        D_H_benefit_sum.append(D_H_benefit_sum_temp)
        D_L_benefit_num.append(D_L_benefit_num_temp)
        D_L_benefit_sum.append(D_L_benefit_sum_temp)
        D_H_not_benefit_sum.append(D_H_not_benefit_sum_temp)
        D_L_compare_to_eq_sum.append(D_L_compare_to_eq_sum_temp)
        D_H_benefit_ave.append(D_H_benefit_ave_temp)
        D_L_benefit_ave.append(D_L_benefit_ave_temp)
        # ---------------------------------------------

    save_path = 'matlab_files/results/don_' + str(M) + '_' + str(mu) + '_' + str(theta) + '_' + str(Lambda) + '_' + str(d)  \
                + '_'+str(c_thre1) + '_' + str(c_thre2) + '_' + str(NPI_change) + '_' + str(NPI_change_tau) + '_' + str(NPI_change_strong) \
                +'_' + input_info['vac_max_rate_choice'] + '_' + input_info['if_wanned_natural_immunity'] \
                + '_' + str(VAS_each) + '_' + str(chi) + '.xlsx'

    df = pd.DataFrame(H_benefit_num)
    df.to_excel(save_path, sheet_name='H_benefit_num')

    df = pd.DataFrame(H_benefit_sum)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='H_benefit_sum')

    df = pd.DataFrame(H_not_benefit_sum)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='H_not_benefit_sum')

    df = pd.DataFrame(L_benefit_num)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='L_benefit_num')

    df = pd.DataFrame(L_benefit_sum)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='L_benefit_sum')

    df = pd.DataFrame(L_compare_to_eq_sum)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='L_compare_to_eq_sum')

    df = pd.DataFrame(L_not_benefit_sum)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='L_not_benefit_sum')

    df = pd.DataFrame(H_benefit_ave)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='H_benefit_ave')

    df = pd.DataFrame(L_benefit_ave)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='L_benefit_ave')

    # --------------------  add deaths info ---------------------------------------
    df = pd.DataFrame(D_H_benefit_num)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='D_H_benefit_num')
    df = pd.DataFrame(D_H_benefit_sum)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='D_H_benefit_sum')
    df = pd.DataFrame(D_H_not_benefit_sum)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='D_H_not_benefit_sum')
    df = pd.DataFrame(D_L_benefit_num)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='D_L_benefit_num')
    df = pd.DataFrame(D_L_benefit_sum)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='D_L_benefit_sum')
    df = pd.DataFrame(D_L_compare_to_eq_sum)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='D_L_compare_to_eq_sum')

    df = pd.DataFrame(D_H_benefit_ave)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='D_H_benefit_ave')

    df = pd.DataFrame(D_L_benefit_ave)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='D_L_benefit_ave')
    # -----------------------------------------------------------------------------


def donation_compare_specific(M=4, mu=5.6e-3, theta=0.26, Lambda=100, chi=0.8,
                              tau=365 * 5, Tau=365 * 0.5, d=15, VAS_each=1, VAS=100,
                              c_thre1=0.3, c_thre2=1.0, NPI_change=0, NPI_change_tau=365 * 1, NPI_change_strong=7):
    # DELTA I (0.1, 2E-5) (0.5, 5E-5) (0.9, 8E-5)
    donation_c = {}
    donation_num = {}
    active_I = {}
    end_time = {}

    model_eq = paras.DiscreteSVEIRD(Tau, VAS_each, 0, VAS_each, M, d, mu, theta, K=[0, 0.5, 0], tau=tau,
                                    c_thre1=c_thre1, c_thre2=c_thre2,NPI_change=NPI_change, NPI_change_tau=NPI_change_tau,
                                    NPI_change_strong=NPI_change_strong, Lambda=Lambda)
    model_eq.trans_SVEIRD()
    active_I['eq' + 'H'] = np.array(model_eq.trans_length(model_eq.active_H)) / popu_sum_H
    active_I['eq' + 'cumH'] = np.array(model_eq.trans_length(model_eq.cum_H)) / popu_sum_H
    active_I['eq' + 'L'] = np.array(model_eq.trans_length(model_eq.active_L)) / popu_sum_L
    active_I['eq' + 'cumL'] = np.array(model_eq.trans_length(model_eq.cum_L)) / popu_sum_L
    end_time['eq'] = [model_eq.pandemic_end_time]

    model_ineq = paras.DiscreteSVEIRD(Tau, VAS, chi, VAS_each, M, d, mu, theta, K=[0, 0.5, 0], tau=tau,
                                 c_thre1=c_thre1, c_thre2=c_thre2,NPI_change=NPI_change, NPI_change_tau=NPI_change_tau,
                                 NPI_change_strong=NPI_change_strong, Lambda=Lambda)
    model_ineq.trans_SVEIRD()
    active_I['ineq' + 'H'] = np.array(model_ineq.trans_length(model_ineq.active_H)) / popu_sum_H
    active_I['ineq' + 'cumH'] = np.array(model_ineq.trans_length(model_ineq.cum_H)) / popu_sum_H
    active_I['ineq' + 'L'] = np.array(model_ineq.trans_length(model_ineq.active_L)) / popu_sum_L
    active_I['ineq' + 'cumL'] = np.array(model_ineq.trans_length(model_ineq.cum_L)) / popu_sum_L
    end_time['ineq'] = [model_ineq.pandemic_end_time]

    for (K1, K2) in [(8e-5, 0.1), (5E-5, 0.5), (2E-5, 0.9)]:
        model_ineq_ad = paras.DiscreteSVEIRD(Tau, VAS, chi, VAS_each, M, d, mu, theta, K=[K1, K2, 1], tau=tau,
                                             c_thre1=c_thre1, c_thre2=c_thre2,NPI_change=NPI_change, NPI_change_tau=NPI_change_tau,
                                             NPI_change_strong=NPI_change_strong, Lambda=Lambda)
        model_ineq_ad.trans_SVEIRD()
        active_I['ad' + str(K1) + str(K2) + 'H'] = np.array(
            model_ineq_ad.trans_length(model_ineq_ad.active_H)) / popu_sum_H
        active_I['ad' + str(K1) + str(K2) + 'cumH'] = np.array(
            model_ineq_ad.trans_length(model_ineq_ad.cum_H)) / popu_sum_H
        active_I['ad' + str(K1) + str(K2) + 'L'] = np.array(
            model_ineq_ad.trans_length(model_ineq_ad.active_L)) / popu_sum_L
        active_I['ad' + str(K1) + str(K2) + 'cumL'] = np.array(
            model_ineq_ad.trans_length(model_ineq_ad.cum_L)) / popu_sum_L
        end_time['ad' + str(K1) + str(K2)] = [model_ineq_ad.pandemic_end_time]
        donation_c['ad' + str(K1) + str(K2)] = model_ineq_ad.trans_length(
            model_ineq_ad.donation_c / sum(country_info['Income'] == 1))
        donation_num['ad' + str(K1) + str(K2)] = model_ineq_ad.trans_length(model_ineq_ad.donation_num)

    save_path = 'matlab_files/results/don_cont' + str(M) + '_' + str(mu) + '_' + str(theta) + '_' + str(Lambda) + '_' + str(d) \
                + '_'+str(c_thre1) + '_'+str(c_thre2) + '_' + str(NPI_change) + '_' + str(NPI_change_tau) + '_' + str(NPI_change_strong) \
                + '_' + input_info['vac_max_rate_choice'] + '_' + input_info['if_wanned_natural_immunity'] \
                + '_' + str(VAS_each) + '_' + str(chi) + '.xlsx'

    df = pd.DataFrame(donation_c)
    df.to_excel(save_path, sheet_name='donation_c')

    df = pd.DataFrame(donation_num)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='donation_num')

    df = pd.DataFrame(active_I)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='active_I')

    df = pd.DataFrame(end_time)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='end_time')


def donation_compare_hop(M=5, mu=1e-3, theta=0.26,Lambda = 100,
                         VAS=100, d=15, chi=0.8, VAS_each=1, tau=365 * 5, Tau=365 * 0.5,
                         c_thre1=0.3, c_thre2=1.0, NPI_change=0, NPI_change_tau=365 * 1, NPI_change_strong=7,
                         K1_all=(8e-5, 6e-5, 4e-5, 2e-5), K2_all=(0.5, 0.6, 0.7, 0.8)):

    cum_I = {}
    end_time = {}
    save_path = 'matlab_files/results/don_hop_' + str(M) + '_' + str(mu) + '_' + str(theta) + '_' + str(Lambda) + '_' + str(d) \
                + '_'+str(c_thre1) + '_' +str(c_thre2) + '_' + str(NPI_change) + '_' + str(NPI_change_tau) + '_' + str(NPI_change_strong) \
                + '_' + input_info['vac_max_rate_choice'] + '_' + input_info['if_wanned_natural_immunity'] \
                + '_' + str(VAS_each) + '_' + str(chi) + '.xlsx'
    for K3 in range(3, 7):
        for i in range(len(K1_all)):
            K1 = K1_all[i]
            K2 = K2_all[i]
            model_ineq_ad = paras.DiscreteSVEIRD(Tau, VAS, chi, VAS_each, M, d, mu, theta, K=[K1, K2, K3], tau=tau,
                                                 c_thre1=c_thre1, c_thre2=c_thre2, NPI_change=NPI_change, NPI_change_tau=NPI_change_tau,
                                                 NPI_change_strong=NPI_change_strong, Lambda=Lambda)
            model_ineq_ad.trans_SVEIRD()
            cum_I['ad' + str(K1) + str(K2) + str(K3) + 'H'] = np.array(
                model_ineq_ad.trans_length1(model_ineq_ad.D_H)) / popu_sum_H
            cum_I['ad' + str(K1) + str(K2) + str(K3) + 'L'] = np.array(
                model_ineq_ad.trans_length(model_ineq_ad.D_L)) / popu_sum_L
            end_time['ad' + str(K1) + str(K2) + str(K3)] = [model_ineq_ad.pandemic_end_time]

    df = pd.DataFrame(cum_I)
    df.to_excel(save_path, sheet_name='cum_I')

    df = pd.DataFrame(end_time)
    with pd.ExcelWriter(save_path, engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='end_time')


#fig2-3 & S2-S5 & S15 & S18-S19
#run_overall_results(M=5, mu=0.0056, theta=0.2, Lambda=500)

# fig 5
#donation_compare_hop(M=5, mu=5.6e-3, theta=0.2, Lambda=500, d=15, K1_all=(8e-5, 6e-5, 4e-5), K2_all=(0.46, 0.6, 0.8))

# fig 4
# donation_compare(M=5, mu=5.6e-3, theta=0.2, Lambda=500)
# donation_compare_specific(M=5, mu=5.6e-3, theta=0.2, Lambda=500)

# # figS6-figS14 & S16-S17
# run_overall_results(M=6, mu=5.6e-3, theta=0.26, Lambda=100, d=15)
# run_overall_results(M=4, mu=5.6e-6, theta=0.4, Lambda=100, d=15)
# run_overall_results(M=7, mu=5.6e-4, theta=0.2, Lambda=10000, d=15)
# run_overall_results(M=3, mu=5.6e-5, theta=0.5, Lambda=10000, d=15)
# run_overall_results(M=5, mu=5.6e-4, theta=0.1, Lambda=1000, d=15)
# run_overall_results(M=4, mu=5.6e-3, theta=0.22, Lambda=100, d=15)
# run_overall_results(M=3, mu=5.6e-5, theta=0.3, Lambda=10000, d=15)
# run_overall_results(M=10, mu=5.6e-3, theta=0.12, Lambda=1000, d=15)
# run_overall_results(M=9, mu=5.6e-5, theta=0.1, Lambda=100, d=15)


# fig S20
# M, mu, theta, Lambda = 5, 5.6e-4, 0.1, 1000
# donation_compare(M=M, mu=mu, theta=theta, Lambda=Lambda)
# donation_compare_specific(M=M, mu=mu, theta=theta, Lambda=Lambda)

# fig S21
# M, mu, theta, Lambda = 4, 5.6e-3, 0.22, 100
# donation_compare(M=M, mu=mu, theta=theta, Lambda=Lambda)
# donation_compare_specific(M=M, mu=mu, theta=theta, Lambda=Lambda)

# fig S22
# M, mu, theta, Lambda = 3, 5.6e-5, 0.3, 10000
# donation_compare(M=M, mu=mu, theta=theta, Lambda=Lambda)
# donation_compare_specific(M=M, mu=mu, theta=theta, Lambda=Lambda)

# fig S23-S24
# for _, (M, mu, theta, Lambda) in enumerate([(10, 5.6e-3, 0.12, 1000),
#                                              (9, 5.6e-5, 0.1, 100)]):
#
#     donation_compare(M=M, mu=mu, theta=theta, Lambda=Lambda)
#     donation_compare_specific(M=M, mu=mu, theta=theta, Lambda=Lambda)


# # natural immunity fig S25 [1]
# # input change to 'Y'
# M, mu, theta, Lambda = 5, 5.6e-3, 0.2, 500
# run_overall_results(M=M, mu=mu, theta=theta, Lambda=Lambda, d=15)

# # figure S26 [1]
# donation_compare_specific(M=M, mu=mu, theta=theta, Lambda=Lambda)
# donation_compare(M=M, mu=mu, theta=theta, Lambda=Lambda)


# input change to 'N'
# # fig S28-S29
# run_overall_results(M=5, mu=5.6e-3, theta=0.2, Lambda=500, d=15, c_thre1=0.2,c_thre2=0.8)
# donation_compare(M=5, mu=5.6e-3, theta=0.2, Lambda=500, d=15, c_thre1=0.2,c_thre2=0.8)
# donation_compare_specific(M=5, mu=5.6e-3, theta=0.2, Lambda=500, d=15, c_thre1=0.2,c_thre2=0.8)
# # fig S30-S31
# run_overall_results(M=5, mu=5.6e-3, theta=0.2, Lambda=500, d=15, c_thre1=0.5,c_thre2=1.2)
# donation_compare(M=5, mu=5.6e-3, theta=0.2, Lambda=500, d=15, c_thre1=0.5,c_thre2=1.2)
# donation_compare_specific(M=5, mu=5.6e-3, theta=0.2, Lambda=500, d=15, c_thre1=0.5,c_thre2=1.2)
# # fig S32
# run_overall_results(M=3, mu=5.6e-3, theta=0.5, Lambda=500, d=15, c_thre1=0.2,c_thre2=1.4)






