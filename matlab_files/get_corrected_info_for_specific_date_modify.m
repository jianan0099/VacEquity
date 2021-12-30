function [corrected_R,corrected_D,corrected_active] = get_corrected_info_for_specific_date_modify(c_date,c_all,T,e_sample_results,p_sample_results,q_MS_sample_results,q_MN_sample_results,P_MS_unt_sample_results,HorL, if_modify_ep)


static_nu = get_trunc_beta(0.9, 0.1^2,0,1);

reported_I = [];
reported_R = [];
reported_D = [];
reported_active = [];

corrected_I = [];
corrected_I_high = [];
corrected_I_low  = [];

corrected_R = [];
corrected_R_low = [];
corrected_R_high = [];

corrected_D = [];
corrected_D_low = [];
corrected_D_high = [];

corrected_active = [];
corrected_active_low = [];
corrected_active_high = [];

for c_index=1:length(c_all)

    c = c_all(c_index);
    disp(c_index)
    disp(c)
    
    info = T((strcmp(T.iso_code,c) & T.date==c_date),:);
    N = info.population;
    N_te = min(info.estimated_tests, N);
    N_te_p = info.totalCases; 
    R_te_p = info.totalRecoveries;
    D_te_p = info.totalDeaths;
    test_prop = N_te/N; 
    
    reported_I = [reported_I, N_te_p];
    reported_R = [reported_R, R_te_p];
    reported_D = [reported_D, D_te_p];
    reported_active = [reported_active, N_te_p-R_te_p-D_te_p];
   
    if R_te_p == 0
        r_de_lower_bound = 1;
    else
        r_de_lower_bound = R_te_p/(N_te_p-D_te_p);
    end
    
    if r_de_lower_bound==1
        r_sample_results = ones(length(e_sample_results),1);
    else
        r_sample_results = get_trunc_beta(0.9+0.1*r_de_lower_bound,(0.1-0.1*r_de_lower_bound)^2,r_de_lower_bound,1); %0.01
    end
   

    N_actual_infected_total = ones(length(e_sample_results),1);
    N_actual_recovered_total = ones(length(e_sample_results),1);
    N_actual_death_total = ones(length(e_sample_results),1);
    N_actual_active_total = ones(length(e_sample_results),1);


    for i=1:length(e_sample_results)
            
        e = e_sample_results(i);
        p = p_sample_results(i);
        q_MS = q_MS_sample_results(i);
        q_MN = q_MN_sample_results(i);
        P_MS_unt = P_MS_unt_sample_results(i);
        r = r_sample_results(i);
        

        if if_modify_ep==1
            if test_prop>1
                e = 1 - power(1-e, test_prop);
                p = 1 - power(1-p, test_prop);
            end
        end
        p = max(p, 1-(N_te_p + (N - N_te) * (min(N_te_p/N_te, 1)) * ((q_MS - q_MN) * P_MS_unt+q_MN))/N);
        
        N_actual_I = N_actual_infected_sample(e,p,q_MS,q_MN,P_MS_unt,N, N_te,N_te_p);

        N_actual_I_te = (N_te_p-(1-p)*N_te)/(e+p-1);


        fp = N_te_p-e*N_actual_I_te;
 
        max_nu = min(1, R_te_p/(r*fp));
        mean_nu = min(max((R_te_p/r-e*N_actual_I_te+D_te_p)/fp, 0),max_nu);
        

        if mean_nu == 0
             nu = static_nu(i);    
        else
            if mean_nu ==1
               nu = 1;
            else
                if R_te_p < 1 || fp < 1
                   nu = 0;  
                else
                    if max_nu == mean_nu
                        nu = max_nu;
                    else
                    [a,b] = find_beta_shape_params(0.9*max_nu+0.1*mean_nu, (0.1*(max_nu-mean_nu)^2));
                    pd= truncate(makedist('Beta',a,b),mean_nu, max_nu);
                    nu = random(pd);
                    end
                end
            end
        end


        frac = N_actual_I/(e*N_actual_I_te);
        N_actual_infected_total(i) = N_actual_I;

        if abs(R_te_p/r - nu * fp)<1
            N_actual_recovered_total(i) = 0;
        else
            N_actual_recovered_total(i) = (R_te_p/r - nu * fp) * frac;
        end
        
        N_actual_death_total(i) = D_te_p*frac;
        N_actual_active_total(i) =N_actual_infected_total(i)-N_actual_recovered_total(i)-N_actual_death_total(i);

    end

    corrected_I = [corrected_I,median(N_actual_infected_total)];
    corrected_I_high = [corrected_I_high,prctile(N_actual_infected_total,97.5)];
    corrected_I_low  = [corrected_I_low,prctile(N_actual_infected_total,2.5)];
    
    corrected_R = [corrected_R,median(N_actual_recovered_total)];
    corrected_R_low = [corrected_R_low, prctile(N_actual_recovered_total,2.5)];
    corrected_R_high = [corrected_R_high, prctile(N_actual_recovered_total,97.5)];

    corrected_D = [corrected_D,median(N_actual_death_total)];
    corrected_D_low = [corrected_D_low, prctile(N_actual_death_total, 2.5)];
    corrected_D_high = [corrected_D_high, prctile(N_actual_death_total, 97.5)];
    
    corrected_active = [corrected_active,median(N_actual_active_total)];
    corrected_active_low = [corrected_active_low, prctile(N_actual_active_total, 2.5)];
    corrected_active_high = [corrected_active_high, prctile(N_actual_active_total, 97.5)];
    disp(median(N_actual_active_total))

end


corrected_active(abs(corrected_active)<1) = 0;
corrected_IFR = corrected_D./corrected_I;

corr_T = table(c_all,reported_I', reported_R', reported_D', reported_active', corrected_R', corrected_R_low', corrected_R_high', corrected_D',  corrected_D_low', corrected_D_high', corrected_IFR', corrected_active',corrected_active_low', corrected_active_high', corrected_I',corrected_I_high',corrected_I_low', ...
'VariableNames',{'iso_code','reported_I', 'reported_R', 'reported_D', 'reported_active',...
'corrected_R','corrected_R_low', 'corrected_R_high', 'corrected_D', 'corrected_D_low', 'corrected_D_high', 'corrected_IFR', 'corrected_active','corrected_active_low', 'corrected_active_high', 'corrected_I','corrected_I_high','corrected_I_low'});
writetable(corr_T,strcat(HorL,datestr(c_date),'.csv'),'Delimiter',' ') 

end

