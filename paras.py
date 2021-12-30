# set model para
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import hyper_paras_setting
from scipy.sparse import diags
import time
import json
import modified_data_collections as mdc

# input ------------------------------------
with open('input.json', 'r') as f:
    input_info = json.load(f)
hyper_paras = hyper_paras_setting.hyper()
# ------------------------------------------

# ------ read modified country info --------------------------------------
country_info = mdc.read_specific_date_info(hyper_paras['t0_date_owid'])
# ------------------------------------------------------------------------

c_index = list(country_info['common']).index('IDN')


def one_hop_further(adj):
    Indenti = np.diag(np.ones(len(adj)))
    adj_temp = adj + np.diag(np.ones(len(adj))) + Indenti
    return (adj_temp.dot(country_info['P']) > 0) - Indenti


def neighbor_based_reallocate(type):
    P1 = country_info['P'] > 0
    P2 = one_hop_further(P1)
    P3 = one_hop_further(P2)
    P4 = one_hop_further(P3)
    if type == 0:
        return country_info['P']
    if type == 1:
        return P1
    if type == 2:
        return P2
    if type == 3:
        return P3
    if type == 4:
        return P4


def neighbor_based_prop(P_):
    LMIC_N = np.diag(country_info['N']).dot(P_)
    LMIC_N = np.transpose(LMIC_N[country_info['Income'] == 0])[country_info['Income'] == 1]
    LMIC_N = LMIC_N / LMIC_N.sum(axis=1)[:, None]
    return LMIC_N


LMIC_N = {}
for type_index in range(5):
    LMIC_N[type_index] = neighbor_based_prop(neighbor_based_reallocate(type_index))


def get_single_part_vaccine_distribution(total_vaccines, allocation_frac, max_need):
    big_omega_temp = total_vaccines * allocation_frac
    big_omega = np.minimum(big_omega_temp, max_need)
    remain_supply = max_need - big_omega
    all_remain_vaccines = total_vaccines - np.sum(big_omega)

    if np.sum(remain_supply) > 0 and all_remain_vaccines > 0:
        further_supply = all_remain_vaccines * remain_supply / np.sum(remain_supply)
        big_omega = np.minimum(big_omega + further_supply, max_need)
    return big_omega


def define_alloc_frac(**kwargs):
    if kwargs['alloc_strategy'] == 1:
        proportional_part = np.array(kwargs['N'])
        allocation_frac = proportional_part / np.sum(proportional_part)
    if kwargs['alloc_strategy'] == 4:
        proportional_part = np.array(kwargs['cum']) \
                            - np.array(kwargs['cum_prime'])
        allocation_frac = proportional_part / np.sum(proportional_part) if np.sum(
            proportional_part) > 0 else np.zeros_like(proportional_part)
    if kwargs['alloc_strategy'] == 6:
        proportional_part = np.array(kwargs['I'])
        allocation_frac = proportional_part / np.sum(proportional_part) if np.sum(
            proportional_part) > 0 else np.zeros_like(proportional_part)

    if kwargs['alloc_strategy'] == 7:
        proportional_part = np.array(kwargs['cum_D']) \
                            - np.array(kwargs['cum_D_prime'])
        allocation_frac = proportional_part / np.sum(proportional_part) if np.sum(
            proportional_part) > 0 else np.zeros_like(proportional_part)
    # ------------------------------------------------------
    return allocation_frac


class VaccineAllocationFuncs:
    """vac distribution"""

    def __init__(self, tau, Tau, VAS, chi, VAS_each):
        # share of global popu that have fully vaccinated at the beginning of the simulation
        # tau: model works in a daily time step, tau is the total time length
        # omega: share of global popu that have fully vaccinated at the end of the first year
        self.varphi_0 = np.sum(country_info['V']) / np.sum(country_info['N'])
        self.vaccine_allocation_strategy = VAS  # vaccine allocation strategy
        self.chi = chi  # vaccine allocation prop in rich countries
        self.vaccine_allocation_strategy_in_each_part = VAS_each
        self.v = math.exp((math.log(0.5) - math.log(self.varphi_0)) / Tau) - 1
        self.varphi_t = []
        for t in range(tau + 1):
            total = self.varphi_0 * pow(1 + self.v, t)
            if total <= 1:
                self.varphi_t.append(total)
            else:
                self.varphi_t.append(
                    self.varphi_t[-1] + self.varphi_0 * pow(1 + self.v, Tau) - self.varphi_0 * pow(1 + self.v, Tau - 1))
        self.delta_varphi = [self.varphi_t[t + 1] - self.varphi_t[t] for t in
                             range(tau)]  # vaccine production at time t=  varphi(t+1)-varphi(t)
        self.sum_popu = np.sum(country_info['N'])

    def vaccine_allocation(self, t, S, remain_vaccines, **kwargs):
        max_need = np.maximum(S - remain_vaccines,
                              np.zeros_like(S, dtype=float))

        global_vaccine_available = min(self.sum_popu * self.delta_varphi[t],
                                       np.sum(max_need))

        if self.vaccine_allocation_strategy == 1 or \
                (self.vaccine_allocation_strategy == 100 and
                 self.vaccine_allocation_strategy_in_each_part == 1):
            allo_frac = define_alloc_frac(alloc_strategy=1, N=kwargs['N'])
        if self.vaccine_allocation_strategy == 4 or \
                (self.vaccine_allocation_strategy == 100 and
                 self.vaccine_allocation_strategy_in_each_part == 4):
            allo_frac = define_alloc_frac(alloc_strategy=4, cum=kwargs['cum'], cum_prime=kwargs['cum_prime'])
        if self.vaccine_allocation_strategy == 6 or \
                (self.vaccine_allocation_strategy == 100 and
                 self.vaccine_allocation_strategy_in_each_part == 6):
            allo_frac = define_alloc_frac(alloc_strategy=6, I=kwargs['I'])

        if self.vaccine_allocation_strategy == 7 or \
                (self.vaccine_allocation_strategy == 100 and
                 self.vaccine_allocation_strategy_in_each_part == 7):
            allo_frac = define_alloc_frac(alloc_strategy=7, cum_D=kwargs['cum_D'], cum_D_prime=kwargs['cum_D_prime'])
        # -------------------------------------------------------------------

        if self.vaccine_allocation_strategy != 100:
            big_omega_results = get_single_part_vaccine_distribution(global_vaccine_available,
                                                                     allocation_frac=allo_frac,
                                                                     max_need=max_need)
        else:
            big_omega_results = np.zeros_like(S, dtype=np.float)
            HC_frac_temp = np.sum(allo_frac[country_info['Income'] == 1])
            HC_frac_total = max(self.chi, HC_frac_temp)
            HC_vaccine_available = global_vaccine_available * HC_frac_total

            HC_frac = allo_frac[country_info['Income'] == 1] / HC_frac_temp if HC_frac_temp > 0 else np.zeros_like(
                allo_frac[country_info['Income'] == 1])
            big_omega_results[country_info['Income'] == 1] = get_single_part_vaccine_distribution(
                HC_vaccine_available, allocation_frac=HC_frac, max_need=max_need[country_info['Income'] == 1])

            LC_vaccine_available = global_vaccine_available * (1 - HC_frac_total)
            LC_frac = allo_frac[country_info['Income'] == 0] / (1 - HC_frac_temp) if (1 - HC_frac_temp) > 0 \
                else np.zeros_like(allo_frac[country_info['Income'] == 0])
            big_omega_results[country_info['Income'] == 0] = get_single_part_vaccine_distribution(
                LC_vaccine_available, allocation_frac=LC_frac, max_need=max_need[country_info['Income'] == 0])
        return big_omega_results

    def compare_vac_model_and_real(self):
        real = np.array(pd.read_excel("raw_data/world_3.1_5.13.xlsx", header=None))[:, 1] / 100
        plt.plot(real, label='real')
        plt.plot(self.varphi_t, label='model')
        plt.legend()
        plt.show()


def get_phi_with_S(big_phi_current, S_current, dv=hyper_paras["dv"]):
    """
    the number of susceptible individuals in each country become fully vaccinated 1
    :param big_phi_current:  number of individuals that can be fully vaccinated with vaccines available at t
    :param S_current: num of susceptible individuals at time t
    :param dv: the time needed for the body to build full immunity against the virus
    :return:
    """
    return np.true_divide(np.minimum(np.minimum(big_phi_current, country_info['Vac_rate_max']), S_current), dv)


class DiscreteSVEIRD:
    """
    for disease transmission
    """

    def __init__(self, Tau, VAS, chi, VAS_each, M, d, mu, theta, K, t_prime=15, tau=365 * 5, gamma=0.00015,
                 c_thre1=0.3, c_thre2=1, NPI_change=0, NPI_change_tau=365 * 1, NPI_change_strong=7, Lambda=5):

        self.tau = tau  # model works in a daily time step, tau is the total time length
        self.VacAll = VaccineAllocationFuncs(tau, Tau, VAS, chi, VAS_each)  # vac related model
        self.t_prime = t_prime  # for calculate trans speed [days]
        self.gamma = gamma

        self.P = country_info['P']
        self.Income = country_info['Income']
        self.N0 = country_info['N']
        self.V0 = country_info['V']
        self.R0 = country_info['R']
        self.D0 = country_info['D']
        self.I0 = country_info['IS']
        self.S0 = self.N0 - self.V0 - self.I0 - self.R0 - self.D0
        self.common = country_info['common']

        self.sum_popu = np.sum(self.N0)
        self.c_num = len(self.N0)
        self.M = M
        self.T = hyper_paras['R0_basic'] * hyper_paras['alpha'] * np.diag(
            np.array([(1 + theta) ** i for i in range(M)]))

        self.U = diags([[1 - mu / (Lambda ** i) for i in range(M - 1)] + [1],
                        [mu / (Lambda ** i) for i in range(M - 1)]], [0, 1]).toarray()
        self.eta = hyper_paras["eta_base"] * np.array([math.exp(-(i / d) ** 2) for i in range(M)])
        self.epsilon = hyper_paras["epsilon_base"] * np.array([math.exp(-(i / d) ** 2) for i in range(M)])
        self.K = K
        self.M_eta = np.diag(1 - self.eta)


        self.F = np.diag(country_info['CFR']).dot(
            np.tile(np.array([(1 + theta) ** i for i in range(M)]) * (1 - self.eta) / (1 - self.eta[0]), (
            len(country_info['CFR']),
            1)))


        self.c_thre1 = c_thre1
        self.c_thre2 = c_thre2

        self.NPI_change = NPI_change
        self.NPI_change_tau = NPI_change_tau
        self.c_strong1 = 1.1 * self.N0 / (self.S0 + self.V0) / NPI_change_strong
        self.c_mild1 = 1.5 * self.N0 / (self.S0 + self.V0) / NPI_change_strong

        # -------------------------------------------------------------------

        self.sigma = hyper_paras["sigma"]
        self.alpha = hyper_paras['alpha']
        self.varepsilon = hyper_paras['varepsilon']

        self.infected_number_IS = []
        self.infected_number_IV = []
        self.infected_number_ES = []
        self.infected_number_EV = []
        self.total_uninfected_number = []
        self.total_death = []
        self.cum = []
        self.cum1 = []
        self.active_H = []
        self.active_L = []

        self.D_H = []  # H
        self.D_L = []  # L
        # -----------------------------------


        self.c_Specific = []
        self.IS_Specific = []
        self.Re_Specific = []
        self.Re_H = []
        self.Re_L = []
        # -------------------------------------


        self.SV_specific = []
        self.VS_specific = []
        # ------------------------------------

        self.cum_all = np.array(self.D0 + self.R0 + self.I0)  # np.zeros(self.c_num)

        self.cum_D_all = np.array(self.D0)
        # ------------------------------------------------------
        self.vac_L = []
        self.vac_H = []
        self.donation_c = []
        self.donation_num = []
        self.cum_H = []
        self.cum_L = []
        self.strain_H = []
        self.strain_L = []

        self.H_to_H = []
        self.L_to_H = []
        # ------------------------------------------


        self.HL_frac_H = []
        self.HL_frac_L = []
        # --------------------------------------------

    def model_mobility_trans(self, X):
        return np.transpose(np.tile(X * self.gamma, (self.c_num, 1))) * self.P

    def reset_allo_c(self, active_case_rate, current_reallo_c, allo_ratio):
        reallocate = np.zeros(self.c_num, dtype=np.int)
        xx = {}
        for i in range(self.c_num):
            if current_reallo_c[i]:
                xx[i] = active_case_rate[i]
        I_sort = np.array(list(xx.values())).argsort()[:math.ceil(sum(current_reallo_c) * allo_ratio)]
        reallocate[np.array(list(xx.keys()))[I_sort]] = 1
        return reallocate

    def cal_Re(self, c_level, S, V, I_sum):
        X = 1 / self.alpha * np.diag(c_level / self.N0).dot(
            (np.tile(S, (self.M, 1)).transpose() + np.tile(V, (self.M, 1)).transpose().dot(self.M_eta))).dot(
            self.T * np.diag(np.diagonal(self.U)))
        #X[I_sum < 1] = 0
        Rt = X.max(axis=1)
        return Rt

    # -----------------------------------------------------------


    def set_c2(self, Re):
        c_current = np.zeros(self.c_num)
        c_current[Re < self.c_thre1] = 1 - 0.4
        c_current[(Re >= self.c_thre1) & (Re < self.c_thre2)] = 1 - 0.4
        c_current[Re >= self.c_thre2] = 1 - 0.8
        return c_current

    @staticmethod
    def E_out_in(E, all_mobility, could_move_nodes):
        E_frac = E / np.transpose(np.tile(could_move_nodes, (E.shape[1], 1)))
        c_out = np.diag(np.sum(all_mobility, axis=1)).dot(E_frac)
        c_in = (all_mobility / np.tile(np.sum(all_mobility, axis=1), (all_mobility.shape[0], 1))).dot(c_out)
        return c_out, c_in

    def E_out_in_cont(self, E_all, all_mobility, could_move_nodes):
        E_frac = E_all / np.transpose(np.tile(could_move_nodes, (E_all.shape[1], 1)))
        c = np.sum(all_mobility[:, country_info['Income'] == 1], axis=1) * np.sum(E_frac, axis=1)
        HIC_exposed_import = np.sum(c[country_info['Income'] == 1])
        LMIC_exposed_import = np.sum(c[country_info['Income'] == 0])
        return HIC_exposed_import, LMIC_exposed_import
    # -----------------

    def trans_SVEIRD(self):
        S = self.S0.copy()
        V = self.V0.copy()
        ES = np.zeros((self.c_num, self.M))
        EV = np.zeros((self.c_num, self.M))
        IS = np.zeros((self.c_num, self.M))
        IS[:, 0] = self.I0.copy()
        IV = np.zeros((self.c_num, self.M))
        I = np.sum(IS + IV, axis=1)
        R = self.R0.copy()
        D = self.D0.copy()
        N = self.N0.copy()
        self.cum_V = [0]
        self.mean_vac_rate = []

        # save time
        T_array = np.tile(np.array([self.T[i][i] for i in range(self.M)]), (self.c_num, 1))
        TV_array = np.tile(np.array([self.T[i][i] * (1 - self.eta[i]) for i in range(self.M)]), (self.c_num, 1))


        F_matrix = self.F  # np.tile(self.F, (self.c_num, 1))
        F_V_matrix = self.F * (
            np.tile(1 - self.epsilon, (self.c_num, 1)))  # np.tile(self.F * (1 - self.epsilon), (self.c_num, 1))


        remain_big_phi = np.zeros(self.c_num)

        cum_prime = np.zeros(self.c_num)
        cum_prime_temp = [self.cum_all]


        cum_D_prime = np.zeros(self.c_num)
        cum_D_prime_temp = [self.cum_D_all]


        for t in range(self.tau):
            if self.VacAll.vaccine_allocation_strategy == 1 or \
                    (self.VacAll.vaccine_allocation_strategy == 100 and
                     self.VacAll.vaccine_allocation_strategy_in_each_part == 1):
                big_omega = self.VacAll.vaccine_allocation(t, S=S, N=N,
                                                           remain_vaccines=remain_big_phi)
            if self.VacAll.vaccine_allocation_strategy == 4 or \
                    (self.VacAll.vaccine_allocation_strategy == 100 and
                     self.VacAll.vaccine_allocation_strategy_in_each_part == 4):
                big_omega = self.VacAll.vaccine_allocation(t, cum=self.cum_all / N, cum_prime=cum_prime / N,
                                                           remain_vaccines=remain_big_phi, S=S)
            if self.VacAll.vaccine_allocation_strategy == 6 or \
                    (self.VacAll.vaccine_allocation_strategy == 100 and
                     self.VacAll.vaccine_allocation_strategy_in_each_part == 6):
                big_omega = self.VacAll.vaccine_allocation(t, I=I / N, remain_vaccines=remain_big_phi, S=S)

            if self.VacAll.vaccine_allocation_strategy == 7 or \
                    (self.VacAll.vaccine_allocation_strategy == 100 and
                     self.VacAll.vaccine_allocation_strategy_in_each_part == 7):
                big_omega = self.VacAll.vaccine_allocation(t, cum_D=self.cum_D_all / N, cum_D_prime=cum_D_prime / N,
                                                           remain_vaccines=remain_big_phi, S=S)


            if self.K[2] == 1:
                # if self.VacAll.vaccine_allocation_strategy == 100:
                # reallocate
                reallocate = np.zeros(self.c_num, dtype=np.int)
                reallocate[(self.Income == 1) & (I < self.K[0] * N)] = 1
                self.donation_c.append(sum(reallocate))
                re_vac = np.sum(big_omega[reallocate == 1] * self.K[1])
                self.donation_num.append(re_vac)
                big_omega[reallocate == 1] = big_omega[reallocate == 1] * (1 - self.K[1])
                big_omega[self.Income == 0] = big_omega[self.Income == 0] + re_vac * define_alloc_frac(
                    alloc_strategy=self.VacAll.vaccine_allocation_strategy_in_each_part,
                    I=I[self.Income == 0] / N[self.Income == 0],
                    cum=self.cum_all[self.Income == 0] / N[self.Income == 0],
                    cum_prime=cum_prime[self.Income == 0] / N[self.Income == 0],
                    N=N[self.Income == 0])

            if self.K[2] >= 2 and self.K[2] < 10:  # K[2]-2 -> LMIC_N
                reallocate = np.zeros(self.c_num, dtype=np.int)
                reallocate[(self.Income == 1) & (I < self.K[0] * N)] = 1
                self.donation_c.append(sum(reallocate))
                re_vac_high = big_omega.copy()
                re_vac_high[reallocate == 0] = 0
                re_vac_high = re_vac_high * self.K[1]
                re_vac = np.sum(big_omega[reallocate == 1] * self.K[1])
                self.donation_num.append(re_vac)

                big_omega[reallocate == 1] = big_omega[reallocate == 1] * (1 - self.K[1])
                big_omega[self.Income == 0] = big_omega[self.Income == 0] + \
                                              np.sum(np.diag(re_vac_high[self.Income == 1]).dot(LMIC_N[self.K[2] - 2]),
                                                     axis=0)

            if self.K[2] >= 10:  # K[2]-10 I mininum
                reallocate = self.reset_allo_c(I / N, (self.Income == 1) & (I < self.K[0] * N), self.K[2] - 10)
                self.donation_c.append(sum(reallocate))
                re_vac = np.sum(big_omega[reallocate == 1] * self.K[1])
                self.donation_num.append(re_vac)
                big_omega[reallocate == 1] = big_omega[reallocate == 1] * (1 - self.K[1])
                big_omega[self.Income == 0] = big_omega[self.Income == 0] + re_vac * define_alloc_frac(
                    alloc_strategy=self.VacAll.vaccine_allocation_strategy_in_each_part,
                    I=I[self.Income == 0] / N[self.Income == 0],
                    cum=self.cum_all[self.Income == 0] / N[self.Income == 0],
                    cum_prime=cum_prime[self.Income == 0] / N[self.Income == 0],
                    N=N[self.Income == 0])


            big_phi = remain_big_phi + big_omega
            self.vac_H.append(sum(V[self.Income == 1]) / sum(country_info['N'][self.Income == 1]))
            self.vac_L.append(sum(V[self.Income == 0]) / sum(country_info['N'][self.Income == 0]))

            phi = get_phi_with_S(big_phi, S)

            self.total_uninfected_number.append(S + V)
            self.total_death.append(D)


            remain_big_phi = big_phi - phi

            I_sum = IS + IV
            if t == 0:
                c_current = country_info['c_init']
                Re = self.cal_Re(c_current,S,V,I_sum)
            else:
                c_current = self.set_c2(Re)
                Re = self.cal_Re(c_current, S, V, I_sum)

            self.Re_H.append(1-np.mean(c_current[self.Income == 1]))
            self.Re_L.append(1-np.mean(c_current[self.Income == 0]))


            self.HL_frac_H.append(np.sum(np.sum(I_sum, axis=1)[self.Income == 1]) / np.sum(I_sum))
            self.HL_frac_L.append(np.sum(np.sum(I_sum, axis=1)[self.Income == 0]) / np.sum(I_sum))
            # -------------------------------------------------

            S_ES_trans = np.diag(S / N * c_current).dot(I_sum * T_array).dot(self.U)
            S_V_trans = phi
            V_S_trans = self.varepsilon * V

            self.cum_V.append(self.cum_V[-1] + np.sum(phi))
            # ------------------------------------
            V_EV_trans = np.diag(V / N * c_current).dot(I_sum * TV_array).dot(self.U)

            ES_IS_trans = self.sigma * ES

            EV_IV_trans = self.sigma * EV

            IS_R_trans = self.alpha * IS * (1 - F_matrix)

            IS_D_trans = self.alpha * IS * F_matrix

            IV_R_trans = self.alpha * IV * (1 - F_V_matrix)

            IV_D_trans = self.alpha * IV * F_V_matrix


            if input_info['if_wanned_natural_immunity'] == 'Y':
                R_S_trans = R / (365 * 2)
            else:
                R_S_trans = np.zeros_like(R)


            could_move_nodes = S + V + np.sum(ES + EV, axis=1) + R

            all_mobility = self.model_mobility_trans(could_move_nodes)
            all_mobility_temp = np.tril(all_mobility, -1) + np.transpose(np.triu(all_mobility, 1))
            all_mobility = all_mobility_temp + all_mobility_temp.T

            ES_out_in = self.E_out_in(ES, all_mobility, could_move_nodes)
            EV_out_in = self.E_out_in(EV, all_mobility, could_move_nodes)
            other_out_in = self.E_out_in(np.array([S, V, R]).transpose(), all_mobility, could_move_nodes)


            dS = - np.sum(S_ES_trans, axis=1) - S_V_trans + V_S_trans - other_out_in[0][:, 0] + other_out_in[1][:,
                                                                                                0] + R_S_trans
            # -----------------------------

            dV = - V_S_trans - np.sum(V_EV_trans, axis=1) + S_V_trans - other_out_in[0][:, 1] + other_out_in[1][:, 1]

            dES = S_ES_trans - ES_IS_trans - ES_out_in[0] + ES_out_in[1]

            dEV = V_EV_trans - EV_IV_trans - EV_out_in[0] + EV_out_in[1]


            self.H_to_H.append(self.E_out_in_cont(ES+EV, all_mobility, could_move_nodes)[0])
            self.L_to_H.append(self.E_out_in_cont(ES+EV, all_mobility, could_move_nodes)[1])
            # --------------------------------------------------------------------

            dIS = ES_IS_trans - IS_R_trans - IS_D_trans

            dIV = EV_IV_trans - IV_R_trans - IV_D_trans


            dR = np.sum(IS_R_trans + IV_R_trans, axis=1) - other_out_in[0][:, 2] + other_out_in[1][:, 2] - R_S_trans
            # -----------------------------

            dD = np.sum(IS_D_trans + IV_D_trans, axis=1)


            S = S + dS
            V = V + dV
            ES = ES + dES
            EV = EV + dEV
            IS = IS + dIS
            IV = IV + dIV
            I = np.sum(IS + IV, axis=1)
            R = R + dR
            D = D + dD
            N = S + V + np.sum(ES + EV, axis=1) + I + R + D

            self.cum_all = self.cum_all + np.sum(ES_IS_trans + EV_IV_trans, axis=1)
            self.cum_H.append(np.sum(self.cum_all[self.Income == 1]))
            self.cum_L.append(np.sum(self.cum_all[self.Income == 0]))
            self.active_H.append(np.sum(I[self.Income == 1]))
            self.active_L.append(np.sum(I[self.Income == 0]))

            self.cum_D_all = np.array(D)
            # ----------

            self.D_H.append(np.sum(D[self.Income == 1]))
            self.D_L.append(np.sum(D[self.Income == 0]))
            # ---------------------------------------------
            self.cum.append(np.sum(ES_IS_trans + EV_IV_trans))
            self.cum1.append(np.sum(ES_IS_trans + EV_IV_trans, axis=0))  # 分析strain分布
            self.strain_H.append(np.sum(ES_IS_trans[self.Income == 1] + EV_IV_trans[self.Income == 1], axis=0))
            self.strain_L.append(np.sum(ES_IS_trans[self.Income == 0] + EV_IV_trans[self.Income == 0], axis=0))

            if t + 1 - self.t_prime < 0:
                cum_prime_temp.append(self.cum_all)
            else:
                cum_prime = cum_prime_temp[(t + 1 - self.t_prime) % self.t_prime]
                cum_prime_temp[(t + 1) % self.t_prime] = self.cum_all

            # -----------
            if t + 1 - self.t_prime < 0:
                cum_D_prime_temp.append(self.cum_D_all)
            else:
                cum_D_prime = cum_D_prime_temp[(t + 1 - self.t_prime) % self.t_prime]
                cum_D_prime_temp[(t + 1) % self.t_prime] = self.cum_D_all
            # --------------------------------------
            IS_temp = IS.copy()
            IS_temp[IS_temp < 1] = 0

            IV_temp = IV.copy()
            IV_temp[IV_temp < 1] = 0

            ES_temp = ES.copy()
            ES_temp[ES_temp < 1] = 0

            EV_temp = EV.copy()
            EV_temp[EV_temp < 1] = 0

            I_current = np.sum(IS_temp + IV_temp + ES_temp + EV_temp, axis=1)

            self.pandemic_end_time = t + 1
            if np.sum(I_current) == 0:
                print(1)
                break

    def trans_length(self, li):
        # list length -> tau
        li = list(li)
        li.extend([0] * (self.tau - len(li)))
        return li

    # ----- trans length for CUM Values -------------
    def trans_length1(self, li):
        # list length -> tau
        li = list(li)
        li.extend([li[-1]] * (self.tau - len(li)))
        return li
    # ------------------------------------------------
