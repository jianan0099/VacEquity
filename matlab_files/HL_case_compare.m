function HL_case_compare(d, c_thre1,c_thre2, NPI_change, NPI_change_tau, NPI_change_strong,...
vac_max_scenario, if_wanned_natural_immunity, VAS_each)
texts = char(97:108);


% --------- set ----------------------------------
strategies = {'eq','ineq0.8','ineq0.9'};
titles = {'Population size based', 'Prevalence based ','Mortality based'};
row_mean = {'Equitable','Inequitable, \chi=0.8','Inequitable, \chi=0.9'};
info_all = {'H_frac', 'H_D_frac','L_frac', 'L_D_frac'};
marker_size = 4;
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
figure('Units', 'centimeters','Position',[0.767291666666667,8.810625,18,5])
legend_all = {'Strain 1','Strain 2','Strain 3','Strain 4', 'Strain 5','Strain 6','Strain 7','Strain 8','Strain 9', 'Strain 10'};
% ------------------------------------------------

figs = 1;
overall_data_path = strcat('results/overall_','5','_', '0.0056', '_', '0.2','_', '500','_',d, '_', c_thre1,'_',c_thre2,'_', NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '.csv');
endtime_data_path = strcat('results/endtime_','5','_', '0.0056', '_', '0.2','_', '500','_', d, '_',c_thre1,'_', c_thre2,'_',NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '.csv');
T_overall = readtable(overall_data_path,'PreserveVariableNames',true);
T_endtime = readtable(endtime_data_path,'PreserveVariableNames',true);

for s=1:3
subplot(1,4,s,'Position', [0.15+(s-1)*0.2, 0.25,0.13,0.5],'Units','normalized')
end_time = T_endtime.(string(strcat(VAS_each,strategies(s))));
col_name1 = string(strcat(VAS_each,strategies(s),'HL_frac_H'));
col_name2 = string(strcat(VAS_each,strategies(s),'HL_frac_L'));
results1 = T_overall.(col_name1) * 100;
results2 = T_overall.(col_name2) * 100;
HL_frac = [results1(1:end_time(1)),results2(1:end_time(1))];
area(HL_frac)
colororder([252/255 141/255 89/255;
            145/255 191/255 219/255])
hold on
set(gca,'FontSize',font_size)
xlim([0 365*5])
xticks(0:365*1:length(results1(:,1)))
xticklabels({'0','1','2','3','4','5'})
ytickformat('percentage')
if s ==3
   legend({'HICs','LMICs'},'Position',[0.734461550968948,0.322683706070287,0.099014110858549,0.303226172954731], 'Units', 'Normalized')
end
xlabel('t [years]','FontSize',font_size)
text(-0.21, 1.16, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
figs = figs +1;

if s==1
text(0.3, 1.18, row_mean(s), 'Units', 'Normalized','FontSize',font_size);
else
text(0.05, 1.18, row_mean(s), 'Units', 'Normalized','FontSize',font_size);
end

if s == 1
   y = ylabel({'Fraction of', 'active cases'});
   set(y,'Units','normalized','FontSize',font_size,'Position',[-0.501315789473684,0.5,0])
end

end

saveas(gcf,strcat('figs/HL_frac.eps'),'epsc')
end