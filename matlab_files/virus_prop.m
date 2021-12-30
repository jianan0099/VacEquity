% S virus compare
function virus_prop(M, theta,d)
path = strcat('results/virus_prop',M, '_',theta,'_',d,'.csv');

virus = readtable(path);
ylabels = {'$T_m$','$\overline{F_m}$','$\eta_m$','$\epsilon_m$'};
font_size = 7;
texts = char(97:100);
          
figure('Units', 'centimeters','Position',[0.767291666666667,8.810625,18,10])
for p=1:4
subplot(2,2,p)
bar(table2array(virus(:,p+1)),'FaceColor',[115/255 115/255 115/255])

ylabel(ylabels(p),'interpreter','latex');

set(gca,'FontSize',font_size)
hold on
xlabel('Strain m');
text(-0.15, 1.15, texts(p), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
end

path = strcat('figs/virus_', M, '_', theta, '.eps');
saveas(gcf,path,'epsc')
end

