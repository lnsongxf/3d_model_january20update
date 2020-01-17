%============================================
%SMOOTHED SHOCKS
startDate=datenum('01-01-1999');
endDate = datenum('01-12-2016');
T=72;%sample length

Date=linspace(startDate,endDate,T);

% load estimation_dataset_quarterly.mat;
% figure('Name','Smoothed vs. Actual output growth','units','normalized','outerposition',[0 0 1 1]);
% 
% hold on;
% plot(Date,oo_.SmoothedVariables.dy_data,'lineWidth',3,'color','black');
% hold on;
% plot(Date,dy_data(2:2+72-1),'--','color','red','lineWidth',3);
% legend('real','model-implied');
%   xlim([startDate endDate])
%   datetick('x','yy','keeplimits');
% 
% set(gca,'FontSize',30);
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'smoothed_dy','-dpdf');


    f=figure;
    f.Name='Counterfactuals';
    f.Units='normalized';
    f.Position=[0 -0.2 1 1];
subplot(5,3,1);
plot(Date,oo_.SmoothedShocks.epsiA,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Productivity');

% subplot(5,2,2);
% plot(Date,oo_.SmoothedShocks.epsiHd,'color','black','lineWidth',5);
%  xlim([startDate endDate])
%  datetick('x','yy','keeplimits');
% title('Housing Depreciation');
% 
% subplot(5,2,3);
% plot(Date,oo_.SmoothedShocks.epsiHk,'color','black','lineWidth',5);
%  xlim([startDate endDate])
%  datetick('x','yy','keeplimits');
% title('Capital Depreciation');


subplot(5,3,2);
plot(Date,oo_.SmoothedShocks.epsiJ,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Housing Preference');

subplot(5,3,3);
plot(Date,oo_.SmoothedShocks.epsiSe,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Corporate Risk');




subplot(5,3,4);
plot(Date,oo_.SmoothedShocks.epsiSB,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Bank Risk');

subplot(5,3,5);
plot(Date,oo_.SmoothedShocks.epsiWe,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Corporate Net Worth');

subplot(5,3,6);
plot(Date,oo_.SmoothedShocks.epsiWb,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Bank Net Worth');

% subplot(7,2,10);
% plot(Date,oo_.SmoothedShocks.epsimarkup_m,'color','black','lineWidth',5);
%  xlim([startDate endDate])
%  datetick('x','yy','keeplimits');
% title('Mortgage Lending Mark-up');

% subplot(7,2,11);
% plot(Date,oo_.SmoothedShocks.epsimarkup_F,'color','black','lineWidth',5);
%  xlim([startDate endDate])
%  datetick('x','yy','keeplimits');
% title('Corporate Lending Mark-up');

subplot(5,3,7);
plot(Date,oo_.SmoothedShocks.epsiEC,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Consumption Preference');


subplot(5,3,8);
plot(Date,oo_.SmoothedShocks.epsiEbF,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Expected Capital Price');

subplot(5,3,9);
plot(Date,oo_.SmoothedShocks.epsimarkup_m,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Mortgage rate markup');

subplot(5,3,10);
plot(Date,oo_.SmoothedShocks.epsimarkup_F,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Corporate rate markup');

subplot(5,3,11);
plot(Date,oo_.SmoothedShocks.epsiECAB,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Bank Capital');
% set(gca,'FontSize',15);

% subplot(5,3,14);
% plot(Date,oo_.SmoothedShocks.epsiBank,'color','black','lineWidth',5);
%  xlim([startDate endDate])
%  datetick('x','yy','keeplimits');
% title('Bank Wide Shock');
% % set(gca,'FontSize',15);

% subplot(5,3,15);
% plot(Date,oo_.SmoothedShocks.meas_dy,'color','black','lineWidth',5);
%  xlim([startDate endDate])
%  datetick('x','yy','keeplimits');
% title('Measurement Error');
% % set(gca,'FontSize',15);



fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'smoothed_shocks','-dpdf');