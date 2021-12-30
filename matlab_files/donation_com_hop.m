function donation_com_hop(M, mu, theta, Lambda, d, c_thre1, c_thre2, NPI_change, NPI_change_tau, NPI_change_strong,...
vac_max_scenario, if_wanned_natural_immunity, vas, chi, K1, K2, if_add_label,max_values,max_values1,texts)
% 
% M = '5';
% mu = '1e-05';
% ini_strain = 'delta';
% pandemic_scenario = 'pessimistic'; %'optimistic'
% vac_max_scenario = 'original';
% if_wanned_natural_immunity = 'N';
% vas = '1';
% chi = '0.8';
% 
% K1 = {'8e-05','6e-05','4e-05', '2e-05'};
% K2 = {'0.5','0.6','0.7', '0.8'};
% if_add_label = 1;
% max_values = [0.09, 0.2];


titles = [];
for kk=1:length(K1)
    specific_K1 = char(K1(kk));
    titles = [titles, strcat('\delta = ',K2(kk), ', I_{thre} = ', specific_K1(1), 'e-5')];
end

path = strcat('results/don_hop_',M,'_', mu, '_', theta,'_', Lambda,'_', d, '_', c_thre1,'_',c_thre2,'_', NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '_', vas, '_', chi, '.xlsx');
T_case =  readtable(path,'Sheet','cum_I','PreserveVariableNames',true);
T_end =  readtable(path,'Sheet', 'end_time','PreserveVariableNames',true);

income = {'H','L'};
figure('Units', 'centimeters','Position',[0.767291666666667,8.810625,18,6.2])
line_width = 1;
font_size = 7;
figs=1;

color_all = [230/255 97/255 1/255;
            253/255 184/255 99/255;
            178/255 171/255 210/255;
            94/255 60/255 153/255;];
line_style = [':', ':', ':','-'];   
marker_style = ['+', 'o', '*','x']; 
marker_size = [10, 6, 6.5, 5];



for row=1:2 
    for col=1:3
        subplot(2,4,4*(row-1)+col,'Position', [0.2+(col-1)*0.2, 0.6 - (row-1) * 0.45, 0.13, 0.3],'Units','normalized')
        K1K2 = strcat('ad', K1(col),K2(col));
        for K3=3:6
            result = table2array(T_case(:,strcat(strcat(K1K2,num2str(K3)), income(row))))*100;
            end_time = table2array(T_end(1,strcat(K1K2,num2str(K3))));
            plot(result(1:end_time(1)),'-', 'LineWidth', line_width, 'MarkerIndices',1:300:length(result),'MarkerSize',marker_size(K3-2),'Marker',marker_style(K3-2), 'color', color_all(K3-2,:))
            hold on
            if end_time(1)<length(result)
                line([end_time(1),end_time(1)],[0,max_values(row)],'linestyle','--','LineWidth', line_width,'color', color_all(K3-2,:));
            end
        end
        ylim([max_values1(row), max_values(row)])
        yticks(max_values1(row):(max_values(row)-max_values1(row))/2:max_values(row))
        

        set(gca,'FontSize',font_size)
        xticks(0:365*1:length(result(:,1)))
        xticklabels({'0','1','2','3','4','5'})
        xlim([365*0 365*5])
        ytickformat('percentage')
        if row==2 
        xlabel('t [years]','FontSize',font_size)
        end
        if row==1
            %legend('1-hop neighbors','2-hop neighbors','3-hop neighbors','4-hop neighbors','Position',[0.111507246376812 + (0.386666666666667-0.08)*(col-1),0.61,0.1662,0.1486])
            %set(gca,'xtick',[])
            title(titles(col),'FontSize',font_size)
        end
        if row==2 && col==1
            text(-0.5,0.5,'Cumulative mortality rate','FontSize',font_size,'Rotation',90,'fontweight','bold','Units','normalized')
        end
       if row==1 && col==3
            text(1.2,0.15,'For HICs','FontSize',font_size,'Rotation',90,'fontweight','bold','Units','normalized')
       end
       if row==2 && col==3
            text(1.2,0.15,'For LMICs','FontSize',font_size,'Rotation',90,'fontweight','bold','Units','normalized')
       end
       if if_add_label==1
        text(-0.2, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
        figs = figs +1;
       end
    end
end

L(1) = plot(nan, nan, '+-','LineWidth', line_width, 'color', color_all(1,:),'MarkerSize',marker_size(1));
L(2) = plot(nan, nan, 'o-','LineWidth', line_width, 'color',color_all(2,:),'MarkerSize',marker_size(2));
L(3) = plot(nan, nan, '*-','LineWidth', line_width, 'color',color_all(3,:),'MarkerSize',marker_size(3));
L(4) = plot(nan, nan, 'x-','LineWidth', line_width, 'color',color_all(4,:),'MarkerSize',marker_size(4));
legend(L, {'1-hop','2-hop','3-hop','4-hop'},'Position',[0.806717974338306,0.349932795698925,0.107352940212278,0.347306989247312],'Units','normalized')
saveas(gcf,strcat('figs/hop_',M,'_', mu, '_', theta,'_', Lambda,'_', d, '_', c_thre1,'_', c_thre2,'_', NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '_', vas, '_', chi, '.eps'),'epsc')
end
