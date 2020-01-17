%============================================
%SMOOTHED SHOCKS
startDate=datenum('01-01-1998');
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
subplot(7,2,1);
plot(Date,oo_.SmoothedShocks.epsiA,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Productivity');

subplot(7,2,2);
plot(Date,oo_.SmoothedShocks.epsiHd,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Housing Depreciation');

subplot(7,2,3);
plot(Date,oo_.SmoothedShocks.epsiHk,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Capital Depreciation');


subplot(7,2,4);
plot(Date,oo_.SmoothedShocks.epsiJ,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Housing Preference');

subplot(7,2,5);
plot(Date,oo_.SmoothedShocks.epsiSe,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Corporate Risk');


subplot(7,2,6);
plot(Date,oo_.SmoothedShocks.epsiSm,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Household Risk');


subplot(7,2,7);
plot(Date,oo_.SmoothedShocks.epsiSB,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Bank Risk');

subplot(7,2,8);
plot(Date,oo_.SmoothedShocks.epsiWe,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Corporate Net Worth');

subplot(7,2,9);
plot(Date,oo_.SmoothedShocks.epsiWb,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Bank Net Worth');

subplot(7,2,10);
plot(Date,oo_.SmoothedShocks.epsimarkup_m,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Mortgage Lending Mark-up');

subplot(7,2,11);
plot(Date,oo_.SmoothedShocks.epsimarkup_F,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Corporate Lending Mark-up');

subplot(7,2,12);
plot(Date,oo_.SmoothedShocks.epsiEC,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Consumption Preference');

subplot(7,2,13);
plot(Date,oo_.SmoothedShocks.epsiECAB,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Bank Capital Depreciation');


% set(gca,'FontSize',15);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'smoothed_shocks','-dpdf');