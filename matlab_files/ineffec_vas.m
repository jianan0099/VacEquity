function ineffec_vas(M, mu, theta, Lambda, d, c_thre1, c_thre2, NPI_change, NPI_change_tau, NPI_change_strong,...
vac_max_scenario, if_wanned_natural_immunity,c_class)
% M = '5';
% mu = '0.001';
% ini_strain = 'delta';
% pandemic_scenario = 'pessimistic'; %'optimistic'
% vac_max_scenario = 'original';
% if_wanned_natural_immunity = 'N';
overall_data_path = strcat('results/overall_',M,'_', mu, '_', theta,'_', Lambda,'_', d, '_', c_thre1,'_',c_thre2,'_', NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '.csv');
endtime_data_path = strcat('results/endtime_',M,'_', mu, '_', theta,'_', Lambda,'_', d, '_', c_thre1,'_', c_thre2,'_', NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '.csv');
T_overall = readtable(overall_data_path,'PreserveVariableNames',true);
T_endtime = readtable(endtime_data_path,'PreserveVariableNames',true);
texts = char(97:112);

% --------- set ----------------------------------
strategies = {'eq','ineq0.8','ineq0.9'};
titles = {'Population size based', 'Prevalence based ','Mortality rate based', 'Incidence based'};
row_mean = {'Equitable','Inequitable, \chi=0.8','Inequitable, \chi=0.9'};
VAS_each = {'1','6','7','4'};
info_all = {'H_frac', 'H_D_frac','L_frac', 'L_D_frac'};
marker_size = 4;
line_width = 1;
font_size = 7;
color_all = [27/255 158/255 119/255;
             217/255 95/255 2/255;
             117/255 112/255 179/255;
              231/255 41/255 138/255;];     
line_style = ['-', '-', '-', '-', '-'];  
makers_style = {'x', '+', '*', 'o'};
marker_sizes = [6.5, 10, 4, 5.5, 4];
figure('Units', 'centimeters','Position',[3.91583333333333,3.30729166666667,18,8.2])
% ------------------------------------------------

figs = 1;
for s=1:3
    subplot(2,4,figs,'Position', [0.15+(figs-1)*0.23, 0.55,0.14,0.25],'Units','normalized')
    for vas=1:4
    col_name = string(strcat(VAS_each(vas),strategies(s),'cum'));
    col_name_time = string(strcat(VAS_each(vas), strategies(s)));
    end_time = T_endtime.(col_name_time);
    result = T_overall.(col_name)*100;
    plot(result,line_style(vas), 'LineWidth', line_width, 'color', color_all(vas,:),'MarkerIndices',1:150:length(result),'MarkerSize',marker_sizes(vas),'Marker',string(makers_style(vas)));
    hold on  
    end
    set(gca,'FontSize',font_size)
    if s==1
        ylabel({'Cumulative incidence','rate woldwide'},'Units','normalized','Position',[-0.65, 0.4, 0],'FontSize',font_size-1)
    end
    xticks(365*0:365*1:365*5)
    xticklabels({'0','1','2','3','4', '5'})
    xlim([365*2 365*5])
    ytickformat('percentage')
    text(-0.2, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
    figs = figs +1;
    
    if s==1
        text(0.2, 1.3, 'Equitable', 'Units', 'Normalized','FontSize',font_size);
    end

    if s==2
        text(0.2, 1.3, {'Inequitable';'\chi=0.8'}, 'Units', 'Normalized','FontSize',font_size);
    end

    if s==3
        text(0.2, 1.3, {'Inequitable'; '\chi=0.9'}, 'Units', 'Normalized','FontSize',font_size);
    end

end

for s=1:3
    subplot(2,4,4+s,'Position', [0.15+(s-1)*0.23, 0.15,0.14,0.25],'Units','normalized')
    for vas=1:4
    col_name = string(strcat(VAS_each(vas),strategies(s),'cum_D'));
    col_name_time = string(strcat(VAS_each(vas), strategies(s)));
    end_time = T_endtime.(col_name_time);
    result = T_overall.(col_name)*100;
    plot(result,line_style(vas), 'LineWidth', line_width, 'color', color_all(vas,:),'MarkerIndices',1:150:length(result),'MarkerSize',marker_sizes(vas),'Marker',string(makers_style(vas)));
    hold on  
    end
    set(gca,'FontSize',font_size)
        if s==1
        ylabel({'Cumulative mortality rate','rate woldwide'},'Units','normalized','Position',[-0.65, 0.4, 0],'FontSize',font_size-1)
        end
    
    xticks(365*0:365*1:365*5)
    xticklabels({'0','1','2','3','4', '5'})
    xlim([365*2 365*5])
    ytickformat('percentage')
    text(-0.2, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
    figs = figs +1;
    xlabel('t [years]','FontSize',font_size)
end

legend(titles,'Position',[0.766402643662163,0.348087384637765,0.227318040603178,0.3934856490701],'NumColumns',1)
saveas(gcf,strcat('figs/vas_compare_',M,'_', mu, '_', theta,'_', Lambda,'_', d, '_', c_thre1,'_',c_thre2,'_', NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '.eps'),'epsc')
end