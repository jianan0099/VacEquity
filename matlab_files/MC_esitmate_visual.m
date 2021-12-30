function MC_esitmate_visual(c_list,median,low,high,reported,y,error,fig)

function y = getmydata2(a)
y = [];
for dd=1:length(a)
s = num2str(a(dd),'%1.4e');
id = find(s=='e');
s = s(1:id-1);
y = [y;str2num(s)];
end
end

%figure('Position',[114,273,1547,492])
marker_size = 4;
line_width = 1;
font_size = 7;

bar(1:length(c_list),median)          
set(gca, 'XTick', 1:length(c_list),'XTickLabel',c_list);
hold on
if error
er = errorbar(1:length(c_list),median,low-median,high-median);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1.5;
hold on 
end
text(1:length(c_list),high,num2str(getmydata2(high),'%.1f'),'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',font_size-3)

plot(reported, 'r*-', 'LineWidth', line_width,'MarkerSize',marker_size)
hold off
set(gca,'FontSize',font_size)

ylabel(y,'FontSize',font_size)
xtickangle(90)
ylim([0 max(high)*1.2])
text(0.7954,0.7737, fig, 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
%saveas(gcf,strcat(strcat('init_compare_',save_name),'.eps'),'epsc')
end

