function donation_full_compare(M, mu, theta, Lambda, d, c_thre1, c_thre2, NPI_change, NPI_change_tau, NPI_change_strong,...
vac_max_scenario, if_wanned_natural_immunity, vas, chi,max_values, if_with_axes)
% M = '5';
% mu = '0.001';
% ini_strain = 'delta';
% pandemic_scenario = 'pessimistic'; %'optimistic'
% vac_max_scenario = 'original';
% if_wanned_natural_immunity = 'N';
% vas = '1';
% chi = '0.8';


D_or_I = 'D_';

path = strcat('results/don_',M,'_', mu, '_', theta,'_', Lambda,'_', d, '_', c_thre1,'_', c_thre2,'_', NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '_', vas, '_', chi, '.xlsx');

T1 = readtable(path, 'Sheet', string(strcat(D_or_I, 'H_benefit_num')));
T2 = readtable(path, 'Sheet', string(strcat(D_or_I, 'H_benefit_ave')));
T3 = readtable(path, 'Sheet', string(strcat(D_or_I, 'L_benefit_num')));
T4 = readtable(path, 'Sheet', string(strcat(D_or_I, 'L_benefit_sum')));

% %T2 = readtable(path, 'Sheet', string(strcat(D_or_I, 'H_benefit_sum')));
T5 = readtable(path, 'Sheet', string(strcat(D_or_I, 'H_not_benefit_sum')));
% %T5 = readtable(path, 'Sheet', string(strcat(D_or_I, 'L_benefit_ave')));
T6 = readtable(path, 'Sheet', string(strcat(D_or_I, 'L_compare_to_eq_sum')));



h = height(T1);
w = width(T1);
T1 = table2array(T1(1:h, 2:w));
T2 = table2array(T2(1:h, 2:w));
T3 = table2array(T3(1:h, 2:w));
T4 = table2array(T4(1:h, 2:w));
T5 = table2array(T5(1:h, 2:w));
T6 = table2array(T6(1:h, 2:w));
T = T1*100;
T(:,:,2) = T2*100;
T(:,:,3) = T3*100;
T(:,:,4) = T4*100;
T(:,:,5) = T5*100;
T(:,:,6) = T6*100;

figure('Units', 'centimeters','Position',[3.91583333333333,3.30729166666667,18,15])
marker_size = 4;
line_width = 1;
font_size = 7;
texts = char(97:108);

x = linspace(0,1, 51);%0,1, 51
y = linspace(0,1e-4, 51);%0,1e-4, 51
[X,Y] = meshgrid(x,y);

figs = 1;
for row=1:2
for col=1:2
    i=(row-1)*2+col;
    ax = subplot(4,2,2*(row-1)+col,'Position', [0.15+(col-1)*0.4, 0.75-(row-1)*0.22, 0.23, 0.16],'Units','normalized');
    t = T(:,:,i);
    t(1,:) = [];
    surf(X,Y,t);
    view(2);
    if row==2 && col==2
        colormap(ax, cbrewer2('YlGn'))
        shading interp
        clb = colorbar;
        t=get(clb,'Limits');
        set(clb,'Ticks',linspace(t(1),t(2),3))
        ylabel(clb, 'r_L','FontSize',font_size)
        clb.Ruler.TickLabelFormat='%.2g%%';
    elseif row==1 && col==2
        colormap(ax, cbrewer2('PiYG'))
        shading interp
        clb = colorbar;
        t=get(clb,'Limits');
        set(clb,'Ticks',linspace(t(1),t(2),3))
        ylabel(clb, 'r_H','FontSize',font_size)
        clb.Ruler.TickLabelFormat='%.2g%%';
    elseif row==1 && col==1
        colormap(ax, cbrewer2('Spectral'))
        shading interp
        clb = colorbar;
        ylabel(clb, {'Fraction of HICs', 'benefiting from donations'},'FontSize',font_size)
        clb.Ruler.TickLabelFormat='%g%%';
    else
        colormap(ax, cbrewer2('Spectral'))
        shading interp
        clb = colorbar;
        ylabel(clb, {'Fraction of LMICs', 'benefiting from donations'},'FontSize',font_size)
        clb.Ruler.TickLabelFormat='%g%%';
    end


    %clb.YRuler.TickLabelFormat = '%.2f';
    
    if row==2
    xlabel('\delta')
    else
    gg = gca;
    gg.XAxis.Visible='off';
    end

    xlim([0.02 1])
    xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    xticklabels({'0','0.2','0.4','0.6','0.8','1'})
    ylim([2e-6 1e-4])
    yticks([0, 2e-5, 4e-5, 6e-5, 8e-5, 10e-5])

    if col==1
        ylabel('I_{thre}')
    else
        set(gca,'ytick',[])
    end
    set(gca,'FontSize',font_size)
    
   
%     if i==1
%         title({'HICs benefiting', 'from donations'},'Position',[0.633853211009174,0.000105,0.00380000006453],'FontSize',font_size,'fontweight','bold');
%     end
%     
%     if row==1 && col==2
%         if strcmp(D_or_I, '')
%             title({'r_H'},'Position',[0.6,0.000105,0.00380000006453],'FontSize',font_size,'fontweight','bold');
%         else
%             title({'r_H'},'Position',[0.6,0.000105,0.00380000006453],'FontSize',font_size,'fontweight','bold');
%         end
% 
%     end
  

%     if row==2 && col==1
%     title({'LMICs benefiting', 'from donations'},'Position',[0.633853211009174,0.000105,0.00380000006453],'FontSize',font_size,'fontweight','bold');
%     end
    
%     if row==2 && col==2
%         if strcmp(D_or_I, '')
%             title({'r_L'},'Position',[0.6,0.000105,0.00380000006453],'FontSize',font_size,'fontweight','bold');
%         else
%              title({'r_L'},'Position',[0.6,0.000105,0.00380000006453],'FontSize',font_size,'fontweight','bold');
%         end
%     %title({'Compared to the inequitable strategy,','for all LMICs'});
%     end
    
%     if i==6
%         if strcmp(D_or_I, '')
%              title({'r_3^L'},'Position',[0.5676,0.000119955957447,0.0038],'FontSize',font_size,'fontweight','bold');
%         else
%             title({'r_3^L'},'Position',[0.733853211009174,0.000110572978724,0.00380000006453],'FontSize',font_size,'fontweight','bold');
%         end
%     %title({'Compared to the equitable strategy,','for all LMICs'});
%     text(1.59,0.000035,1e-15,'For LMICs','FontSize',font_size,'Rotation',90,'fontweight','bold')
%     end

        text(-0.03486238539643,1.25, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
        figs = figs +1;
end
end

path1 = strcat('results/don_cont',M,'_', mu, '_', theta,'_', Lambda,'_', d, '_', c_thre1,'_', c_thre2,'_', NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '_', vas, '_', chi, '.xlsx');
T_donation_c =  readtable(path1,'Sheet','donation_c');
T_donation_num =  readtable(path1,'Sheet','donation_num');
T_case =  readtable(path1,'Sheet','active_I');
T_end =  readtable(path1,'Sheet','end_time');
K1 = {'ad8e_05','ad5e_05','ad2e_05'};
K2 = {'0_1','0_5','0_9'};


color_all = [212/255 185/255 218/255;
            231/255 41/255 138/255;
            142/255 1/255 82/255;
              33/255 102/255 172/255;
              252/255 78/255 42/255;];
line_style = ['-', '-', '-', '-', '-'];   
makers_style = {'x', '+', 's', 'o', ''};
marker_sizes = [15, 12, 8, 11, 20];


subplot(4,3,3*(3-1)+1,'Position', [0.2, 0.3, 0.18, 0.15],'Units','normalized')
for i=1:3
    result = table2array(T_donation_c(:,strcat(K1(i),K2(i))))*100;
    end_time = table2array(T_end(1,strcat(K1(i),K2(i))));
    plot(result(1:end_time(1)),'LineWidth', line_width, 'color', color_all(i,:));
    hold on
    set(gca,'FontSize',font_size)
    xlim([0 365*5])
    xticks(0:365*1:length(result(:,1)))
    xticklabels({'0','1','2','3','4','5'})
end
text(-0.2, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
figs = figs +1;

ylabel({'Fraction of HICs', 'donating vaccines'},'FontSize',font_size,'fontweight','bold');
%yticks(0:0.5:1)
%yticklabels({'0','0.5','1'})
ytickformat('percentage')

subplot(4,3,3*(3-1)+2,'Position', [0.2+(2-1)*0.3, 0.3, 0.18, 0.15],'Units','normalized')
for i=1:3
    result = table2array(T_donation_num(:,strcat(K1(i),K2(i))));
    end_time = table2array(T_end(1,strcat(K1(i),K2(i))));
    plot(result(1:end_time(1)),'LineWidth', line_width, 'color', color_all(i,:));
    hold on
    set(gca,'FontSize',font_size)
    xlim([0 365*5])
    xticks(0:365*1:length(result(:,1)))
    xticklabels({'0','1','2','3','4','5'})
%     if end_time(1)<length(result)
%     line([end_time(1),end_time(1)],[0,4.5e7],'linestyle','--','LineWidth', line_width, 'color', color_all(i,:));
%     end
end

text(-0.2, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
figs = figs +1;

ylabel({'Number of' 'donated vaccines'},'FontSize',font_size,'fontweight','bold','Units', 'Normalized','Position', [-0.162666663090388,0.500000476837158,0]);
%ytickformat('percentage')

subplot(4,3,3*(4-1)+1,'Position', [0.2, 0.08, 0.18, 0.15],'Units','normalized')
for i=1:5
    if i<=3
    result = table2array(T_case(:,strcat(strcat(K1(i),K2(i)),'H')))*100;
    end_time = table2array(T_end(1,strcat(strcat(K1(i),K2(i)))));
    else if i==4
            result = table2array(T_case(:,'eqH'))*100;
            end_time = table2array(T_end(1,'eq'));
        else
            result = table2array(T_case(:,'ineqH'))*100;
            end_time = table2array(T_end(1,'ineq'));
        end
    end
    if i==4
        plot(result(1:end_time(1)),line_style(i), 'LineWidth', line_width, 'color', color_all(i,:),'MarkerIndices',1:50:length(result),'MarkerSize',marker_size,'Marker',string(makers_style(i)))
    else
    plot(result(1:end_time(1)),'LineWidth', line_width, 'color', color_all(i,:));
    end
    hold on
    
    if end_time(1)<length(result)
    line([end_time(1),end_time(1)],[0,max_values(1)],'linestyle','--','LineWidth', line_width,'color', color_all(i,:));
    end
    set(gca,'FontSize',font_size)
    xlim([0 365*5])
    xticks(0:365*1:length(result(:,1)))
    xticklabels({'0','1','2','3','4','5'})
end
xlabel('t [years]','FontSize',font_size)
ylabel({'Prevalence', 'in HICs'},'FontSize',font_size,'fontweight','bold');
ytickformat('percentage')
text(-0.2, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
figs = figs +1;

if if_with_axes==1
if end_time>365*2
axes('Position',[0.2136,0.1018,0.1231,0.0638])
% result = table2array(T_case(:,strcat(strcat(K1(1),K2(1)),'H')))*100;
% end_time = table2array(T_end(1,strcat(strcat(K1(1),K2(1)))));
% plot(result(1:end_time(1)),'LineWidth', line_width, 'color', color_all(1,:));
% hold on
% result = table2array(T_case(:,strcat(strcat(K1(3),K2(3)),'H')))*100;
% end_time = table2array(T_end(1,strcat(strcat(K1(3),K2(3)))));
% plot(result(1:end_time(1)),'LineWidth', line_width, 'color', color_all(3,:));
% 
% result = table2array(T_case(:,'ineqH'))*100;
% end_time = table2array(T_end(1,'ineq'));
% plot(result(1:end_time(1)),'LineWidth', line_width, 'color', color_all(5,:));

for i=1:5
    if i<=3
    result = table2array(T_case(:,strcat(strcat(K1(i),K2(i)),'H')))*100;
    end_time = table2array(T_end(1,strcat(strcat(K1(i),K2(i)))));
    else if i==4
            result = table2array(T_case(:,'eqH'))*100;
            end_time = table2array(T_end(1,'eq'));
        else
            result = table2array(T_case(:,'ineqH'))*100;
            end_time = table2array(T_end(1,'ineq'));
        end
    end
    if i==4
        plot(result(1:end_time(1)),line_style(i), 'LineWidth', line_width, 'color', color_all(i,:),'MarkerIndices',1:50:length(result),'MarkerSize',3,'Marker',string(makers_style(i)))
    else
    plot(result(1:end_time(1)),'LineWidth', line_width, 'color', color_all(i,:));
    end
    hold on
end
xlim([365*2 365*5])
xticks(365*2:365*1:length(result(:,1)))
xticklabels({'2','3','4','5'})
ytickformat('percentage')
set(gca,'FontSize',font_size-3)
end
end
l_all = [];
subplot(4,3,3*(4-1)+2,'Position', [0.2+(2-1)*0.3, 0.08, 0.18, 0.15],'Units','normalized')
for i=1:5
    if i<=3
    result = table2array(T_case(:,strcat(strcat(K1(i),K2(i)),'L')))*100;
    end_time = table2array(T_end(1,strcat(strcat(K1(i),K2(i)))));
    else if i==4
            result = table2array(T_case(:,'eqL'))*100;
            end_time = table2array(T_end(1,'eq'));
        else
            result = table2array(T_case(:,'ineqL'))*100;
            end_time = table2array(T_end(1,'ineq'));
        end
    end
    if i==4
        l = plot(result(1:end_time(1)),line_style(i), 'LineWidth', line_width, 'color', color_all(i,:),'MarkerIndices',1:50:length(result),'MarkerSize',marker_size,'Marker','o');
    else
        l =  plot(result(1:end_time(1)),'LineWidth', line_width, 'color', color_all(i,:));
    end
    hold on
    l_all = [l_all; l];
    if end_time(1)<length(result)
    line([end_time(1),end_time(1)],[0,max_values(2)],'linestyle','--','LineWidth', line_width,'color', color_all(i,:));
    end
    set(gca,'FontSize',font_size)
    xlim([0 365*5])
    xticks(0:365*1:length(result(:,1)))
    xticklabels({'0','1','2','3','4','5'})
end
xlabel('t [years]','FontSize',font_size)
ylabel({'Prevalence','in LMICs'},'FontSize',font_size,'fontweight','bold');
ytickformat('percentage')
text(-0.2, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
figs = figs +1;
for i=1:5
    if i<=3
    end_time = table2array(T_end(1,strcat(strcat(K1(i),K2(i)))));
    else
        if i==4
            end_time = table2array(T_end(1,'eq'));
        else
            end_time = table2array(T_end(1,'ineq'));
        end
    end

end
legend(l_all, '\delta=0.1, I_{thre}=8e-5','\delta=0.5, I_{thre}=5e-5','\delta=0.9, I_{thre}=2e-5','Equitable','Inequitable \chi=0.8','Position',[0.734076693564229,0.13215859030837,0.174746835847536,0.226872246696035])
saveas(gcf,strcat('figs/don_',M,'_', mu, '_', theta,'_', Lambda,'_', c_thre1,'_', c_thre2,'_',NPI_change, '_',NPI_change_tau, '_',NPI_change_strong,...
'_', vac_max_scenario, '_', if_wanned_natural_immunity, '_', vas, '_', chi, '.eps'),'epsc')
end


