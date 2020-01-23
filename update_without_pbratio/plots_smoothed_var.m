%============================================
%SMOOTHED SHOCKS
startDate=datenum('01-01-1999');
endDate = datenum('01-12-2016');
T=72;%sample length

Date=linspace(startDate,endDate,T);



    f=figure;
    f.Name='Counterfactuals';
    f.Units='normalized';
    f.Position=[0 -0.2 1 1];
subplot(5,3,1);
plot(Date,oo_.SmoothedShocks.epsiA,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Productivity');



subplot(5,3,2);
plot(Date,oo_.SmoothedShocks.epsiJ,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Housing Preference');



subplot(5,3,3);
plot(Date,oo_.SmoothedShocks.epsiH,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Housing Price');

subplot(5,3,4);
plot(Date,oo_.SmoothedShocks.epsiSe,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Entrepreneur Risk');

subplot(5,3,5);
plot(Date,oo_.SmoothedShocks.epsiSB,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Bank Risk');

subplot(5,3,6);
plot(Date,oo_.SmoothedShocks.epsiWe,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Entrepreneur Net Worth');

subplot(5,3,7);
plot(Date,oo_.SmoothedShocks.epsimarkup_m,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Mortgage Lending Mark-up');

subplot(5,3,8);
plot(Date,oo_.SmoothedShocks.epsimarkup_F,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Corporate Lending Mark-up');

subplot(5,3,9);
plot(Date,oo_.SmoothedShocks.epsiEC,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Consumption Preference');

subplot(5,3,10);
plot(Date,oo_.SmoothedShocks.epsiECAB,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('EBank Capital');

subplot(5,3,11);
plot(Date,oo_.SmoothedShocks.epsiLTVH,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('HH LTV');

subplot(5,3,12);
plot(Date,oo_.SmoothedShocks.epsiEbF,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Expected Capital Price');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'smoothed_shocks','-dpdf');