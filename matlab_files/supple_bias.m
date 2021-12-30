% modified data compare
figure('Units', 'centimeters','Position',[0.767291666666667,8.810625,18,8])
font_size = 7;
model_start_time = '2021-06-15';
T = readtable(strcat(datestr(model_start_time),'.csv'));
c_all = T.iso_code;

reported_I = T.reported_I;
corrected_I = T.corrected_I;
corrected_I_low = T.corrected_I_low;
corrected_I_high = T.corrected_I_high;

reported_D = T.reported_D;
corrected_D = T.corrected_D;
corrected_D_low = T.corrected_D_low;
corrected_D_high = T.corrected_D_high;

reported_R = T.reported_R;
corrected_R = T.corrected_R;
corrected_R_low = T.corrected_R_low;
corrected_R_high = T.corrected_R_high;

reported_active = T.reported_active;
corrected_active = T.corrected_active;
corrected_active_low = T.corrected_active_low;
corrected_active_high = T.corrected_active_high;

[max_vals, max_index] = maxk(corrected_I,10);
c_list = c_all(max_index);

subplot(2,2,1)
MC_esitmate_visual(c_list,corrected_I(max_index),...
                          corrected_I_low(max_index),...
                          corrected_I_high(max_index),...
                          reported_I(max_index),...
                          'Cumulative infections',1,'a')
title('Estimated numbers ({\color[rgb]{0.00,0.45,0.74}bars}) vs. reported numbers ({\color{red}*})','FontSize',font_size,'Units', 'normalized','Position',[1.110039472579956,1.139150619506836,0])

subplot(2,2,2)
MC_esitmate_visual(c_list,corrected_R(max_index),...
                          corrected_R_low(max_index),...
                          corrected_R_high(max_index),...
                          reported_R(max_index),...
                          'Cumulative recoveries',1,'b')   

subplot(2,2,3)
MC_esitmate_visual(c_list,corrected_D(max_index),...
                          corrected_D_low(max_index),...
                          corrected_D_high(max_index),...
                          reported_D(max_index),...
                          'Cumulative deaths',1,'c')
                      
subplot(2,2,4)     
MC_esitmate_visual(c_list,corrected_active(max_index),...
                          corrected_active_low(max_index),...
                          corrected_active_high(max_index),...
                          reported_active(max_index),...
                          'Active cases',1,'d')     

saveas(gcf,'figs/init_compare.eps','epsc')
