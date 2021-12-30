
function overall_fraction_full(M, mu, theta, Lambda, d, c_thre1, c_thre2,NPI_change, NPI_change_tau, NPI_change_strong,...
vac_max_scenario, if_wanned_natural_immunity, four_limits, four_limits1, yticks_final,time_length)

overall_data_path = strcat('results/overall_',M,'_', mu, '_', theta,'_', Lambda,'_', d, '_', c_thre1,'_', c_thre2,'_', NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '.csv');
endtime_data_path = strcat('results/endtime_',M,'_', mu, '_', theta,'_', Lambda,'_', d, '_', c_thre1,'_', c_thre2,'_', NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '.csv');
T_overall = readtable(overall_data_path,'PreserveVariableNames',true);
T_endtime = readtable(endtime_data_path,'PreserveVariableNames',true);

% --------- set ----------------------------------
strategies = {'eq','ineq0.7','ineq0.8','ineq0.9'};
titles = {'Population size based', 'Prevalence based ','Incidence based','Mortality rate based'};
VAS_each = {'1','6','4','7'};
info_all = {'H_frac', 'H_D_frac','L_frac', 'L_D_frac'};
marker_size = 4;
line_width = 1;
font_size = 7;
color_all = [33/255 102/255 172/255;
              254/255 178/255 76/255;
              252/255 78/255 42/255;
              189/255 0/255 38/255];
line_style = ['-', '-', '-', '-', '-'];   
figure('Units', 'centimeters','Position',[3.91583333333333,3.30729166666667,18,12])

max_values = zeros(20,1);
[max_values(5*(1-1)+1), max_values(5*(1-1)+2), max_values(5*(1-1)+3),  max_values(5*(1-1)+4)] = deal(four_limits(1));
[max_values(5*(2-1)+1), max_values(5*(2-1)+2), max_values(5*(2-1)+3),  max_values(5*(2-1)+4)] = deal(four_limits(2));
[max_values(5*(3-1)+1), max_values(5*(3-1)+2), max_values(5*(3-1)+3),  max_values(5*(3-1)+4)] = deal(four_limits(3));
[max_values(5*(4-1)+1), max_values(5*(4-1)+2), max_values(5*(4-1)+3),  max_values(5*(4-1)+4)] = deal(four_limits(4));

texts = char(97:120);

% ------------------------------------------------

figs = 1;
for row=1:4
    info = info_all(row); % what kind of data
    for col=1:4
        subplot(4,5,5*(row-1)+col,'Position', [0.1+(col-1)*0.19, 0.78-(row-1)*0.22, 0.12, 0.13422],'Units','normalized')
        h_all = [];
        vas = VAS_each(col); % vaccine distribution method
        for i=[1,3,4]
        strategy = strategies(i); % eq or ineq
        col_name_overall = string(strcat(vas, strategy, info));
        col_name_time = string(strcat(vas, strategy));
        end_time = T_endtime.(col_name_time);
        result = T_overall.(col_name_overall) * 100;
        
        if strcmp(strategy,'eq')
        h = plot(result(1:end_time(1)),line_style(i), 'LineWidth', line_width, 'color', color_all(i,:),'MarkerIndices',1:80:length(result),'MarkerSize',marker_size,'Marker','o');
        else
        h = plot(result(1:end_time(1)),line_style(i), 'LineWidth', line_width, 'color', color_all(i,:));
        end
        h_all = [h_all, h];
        
        hold on
        
        if end_time(1)<length(result)
            line([end_time(1),end_time(1)],[0,max_values(5*(row-1)+col)],'linestyle','--','LineWidth', line_width, 'color', color_all(i,:));
        end
        set(gca,'FontSize',font_size)
        if time_length==5
        xlim([0 365*5])
        xticks(0:365*1:length(result(:,1)))
        xticklabels({'0','1','2','3','4','5'})
        else
        xlim([0 365*10])
        xticks(0:365*2:length(result(:,1)))
        xticklabels({'0','2','4','6','8','10'})
        end
        
        
        xlabel('t [years]','FontSize',font_size)
        ytickformat('percentage')
        if row==1
            title(titles(col),'Position',[0.5,1.3,0],'Units','normalized');
        end
        if row==4
            yticks(yticks_final)
        end
        if col==1
        if row == 1
            ylabel({'Prevalence', 'in HICs'},'FontSize',font_size,'Position',[-0.5, 0.5, 0],'Units','normalized')
        end
        if row ==3
            ylabel({'Prevalence', 'in LMICs'},'FontSize',font_size,'Position',[-0.5, 0.5, 0],'Units','normalized')
        end
        if row == 2
            ylabel({'Cumulative mortality','rate in HICs'},'FontSize',font_size,'Position',[-0.5, 0.5, 0],'Units','normalized')
        end
        if row ==4
            ylabel({'Cumulative mortality','rate in LMICs'},'FontSize',font_size,'Position',[-0.5, 0.5, 0],'Units','normalized')
        end
        end
        end
        text(-0.2, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
        figs = figs +1;
        ylim([four_limits1(row) four_limits(row)]);
        yticks(four_limits1(row):(four_limits(row)-four_limits1(row))/2:four_limits(row))
        
    end
end

legend(h_all,'Equitable','Inequitable, \chi=0.8','Inequitable, \chi=0.9','Position',[0.8379,0.3903,0.1421,0.2455])
saveas(gcf,strcat('figs/overall_full_',M,'_', mu, '_', theta,'_', Lambda,'_', c_thre1,'_', c_thre2,'_',NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '.eps'),'epsc')
end

