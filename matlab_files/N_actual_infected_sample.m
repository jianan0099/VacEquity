function N_actual_I = N_actual_infected_sample(e,p,q_MS,q_MN,P_MS_unt,N, N_te,N_te_p)
N_unt = N - N_te;
P_positive_te = min(N_te_p/N_te, 1);
N_actual_I = (N_te_p + N_unt * P_positive_te * ((q_MS - q_MN) * P_MS_unt+q_MN) - (1-p)*N)/(e+p-1);
if abs(N_actual_I)<1
    N_actual_I = 0;
end
end

