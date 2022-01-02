
function strain_compare(M, mu, theta, Lambda, d, c_thre1,c_thre2, NPI_change, NPI_change_tau, NPI_change_strong,...
vac_max_scenario, if_wanned_natural_immunity, ylim1, c_class, VAS_each)
% % strain fraction
% M = '5';
% mu = '0.001';
% ini_strain = 'delta';
% pandemic_scenario = 'pessimistic'; %'optimistic'
% vac_max_scenario = 'original';
% if_wanned_natural_immunity = 'N';
% ylim1 = [0, 0.25];%total [0, 0.15];
% ylim2 = [0, 0.018]; %total [0, 0.008];
% c_class = '';
% VAS_each = {'1'};
% ylabel_add = ' worldwide';%' in LMICs';
overall_data_path = strcat('results/overall_',M,'_', mu, '_', theta,'_', Lambda,'_',d, '_', c_thre1,'_',c_thre2,'_', NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '.csv');
endtime_data_path = strcat('results/endtime_',M,'_', mu, '_', theta,'_', Lambda,'_', d, '_',c_thre1,'_', c_thre2,'_',NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '.csv');
T_overall = readtable(overall_data_path,'PreserveVariableNames',true);
T_endtime = readtable(endtime_data_path,'PreserveVariableNames',true);
texts = char(97:108);

% --------- set ----------------------------------
popu_all = 7660561164;
strategies = {'eq','ineq0.8','ineq0.9'};
titles = {'Population size based', 'Prevalence based ','Mortality based'};
row_mean = {'Equitable','Inequitable, \chi=0.8','Inequitable, \chi=0.9'};
info_all = {'H_frac', 'H_D_frac','L_frac', 'L_D_frac'};
marker_size = 8;
line_width = 1;
font_size = 7;
color_all = [141/255 211/255 199/255;
             190/255 186/255 218/255;
             251/255 128/255 114/255;
             128/255 177/255 211/255;
             253/255 180/255 98/255;         
             179/255 222/255 105/255;
             188/255 128/255 189/255;
             204/255 235/255 197/255;
             217/255 217/255 217/255;
             255/255 237/255 111/255;
             252/255 205/255 229/255; 
             255/255 255/255 179/255;];     
line_style = ['-', '-', '-', '-', '-','-', '-', '-', '-', '-','-', '-'];  
makers_style = {'s', '+', '*', 'o', '.','s', '+', '*', 'o', '.','s', '+'};
marker_sizes = [10, 18, 8, 11, 20,10, 18, 8, 11, 20,10, 18];
figure('Units', 'centimeters','Position',[3.91583333333333,3.30729166666667,18,8.2])
legend_all = {'Strain 1','Strain 2','Strain 3','Strain 4', 'Strain 5','Strain 6','Strain 7','Strain 8','Strain 9', 'Strain 10'};
% ------------------------------------------------

figs = 1;
for s=1:3
subplot(2,4,figs,'Position', [0.15+(figs-1)*0.23, 0.61,0.16,0.3],'Units','normalized')
single_results = [];
for strain=1:str2num(M)
   col_name = strcat(VAS_each(1),strategies(s),string(strain),c_class);
   col_name_time = string(strcat(VAS_each(1), strategies(s)));
   end_time = T_endtime.(col_name_time);
   results = T_overall.(col_name);
   single_results = [single_results,results(1:end_time(1))];           
end
sum_I = sum(single_results,2);
frac_strain =single_results./sum_I * 100;
frac_strain = frac_strain(1:end,:);
area(frac_strain)
colororder(color_all)
hold on
set(gca,'FontSize',font_size)
xlim([0 365*5])
xticks(0:365*1:length(results(:,1)))
xticklabels({'0','1','2','3','4','5'})
ytickformat('percentage')

    if s == 1
       y = ylabel({'Fraction of new cases'});
       set(y,'Units','normalized','FontSize',font_size-1,'Position',[-0.548346126185028,0.5,0],'FontWeight','bold')
    end
    
if s ==3
   legend(legend_all(1:str2num(M)),'Position',[0.818164906220112,0.356453183520599,0.113659021059276,0.3319], 'Units', 'Normalized')
end
text(-0.21, 1.16, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
figs = figs +1;

if s==1
text(0.3, 1.18, row_mean(s), 'Units', 'Normalized','FontSize',font_size);
else
text(0.05, 1.18, row_mean(s), 'Units', 'Normalized','FontSize',font_size);
end


end


for s=1:3
   subplot(2,4,4+s,'Position', [0.15+(s-1)*0.23, 0.71-0.55,0.16,0.3],'Units','normalized')
   for strain=1:str2num(M)
       col_name = strcat(VAS_each(1),strategies(s),string(strain),c_class);
       col_name_time = string(strcat(VAS_each(1), strategies(s)));
       end_time = T_endtime.(col_name_time);
       result = T_overall.(col_name)*100;
       plot(result(1:end_time(1)),line_style(strain), 'LineWidth', line_width, 'color', color_all(strain,:));
       hold on      
       set(gca,'FontSize',font_size)
        xlim([0 365*5])
        xticks(0:365*1:length(result(:,1)))
        xticklabels({'0','1','2','3','4','5'})
        ytickformat('percentage')
   end
    text(-0.21, 1.16, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
    figs = figs +1;
    xlabel('t [years]','FontSize',font_size)
    
    if s == 1
        if strcmp(c_class,'')
            y = ylabel({'Ratio between the'; 'number of new case and';'the world population'});
        elseif c_class == 'L'
            y = ylabel({'Ratio between the'; 'number of new case and';'the population of LMICs'});
        else
            y = ylabel({'Ratio between the'; 'number of new case and';'the population of HICs'});
        end
       set(y,'Units','normalized','FontSize',font_size-1,'Position',[-0.4102,0.5,0],'FontWeight','bold')
    end
   
    ylim(ylim1)
    
end

h = axes('Parent', gcf, 'Position', [0.20660315356237,0.237827715355805,0.079649957301452,0.158724719101124]);
hold(h); % <--- add this line

for strain=1:str2num(M)
   col_name = strcat(VAS_each(1),strategies(1),string(strain),c_class);
   col_name_time = string(strcat(VAS_each(1), strategies(s)));
   end_time = T_endtime.(col_name_time);
   result = T_overall.(col_name)*100;
   plot(result(1:end_time(1)),line_style(strain), 'LineWidth', line_width, 'color', color_all(strain,:));
   hold on      
   set(gca,'FontSize',font_size-1)
    xlim([0 365*1])
    xticks([0 365*0.5 365*1])
    xticklabels({'0','0.5','1'})
    ytickformat('percentage')
end



saveas(gcf,strcat('figs/s_frac_',M,'_', mu, '_', theta,'_', Lambda,'_',  c_thre1,'_',c_thre2,'_', NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '_',  string(VAS_each(1)), '_', c_class,'.eps'),'epsc')
end