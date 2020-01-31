%dy_data, dq_H_data,int_rate_HH_data,int_rate_business_data,bank_rate_data,dbe_data,dbm_data,phib;// , dq_H_data,d_b_to_Y_data,int_rate_HH_data,bank_rate_data,bsp_H_data ,d_b_to_Y_data; 

%run after LTV1.mod
startDate=datenum('01-01-1999');
endDate = datenum('01-12-2016');
T=75;%sample length

Date=linspace(startDate,endDate,T);

figure(1);
% set(gca,'FontSize',30);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'decomp_dy','-dpdf');

figure(2);
% set(gca,'FontSize',30);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'decomp_dqh','-dpdf');


figure(3);
% set(gca,'FontSize',30);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'decomp_int_rate_HH','-dpdf');


figure(4);
% set(gca,'FontSize',30);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'decomp_int_rate_business','-dpdf');

figure(5);
% set(gca,'FontSize',30);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'decomp_bank_rate','-dpdf');


figure(6);
% set(gca,'FontSize',30);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'decomp_dbe','-dpdf');


figure(7);
% set(gca,'FontSize',30);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'decomp_dbm','-dpdf');



figure(8);
% set(gca,'FontSize',30);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'decomp_dinve','-dpdf');

figure(9);
% set(gca,'FontSize',30);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'decomp_dw','-dpdf');

figure(10);
% set(gca,'FontSize',30);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'decomp_dc','-dpdf');



