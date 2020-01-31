startDate=datenum('01-01-1999');
endDate = datenum('01-12-2016');
T=72;%sample length

Date=linspace(startDate,endDate,T);


    f=figure;
    f.Name='Some Smoothed Variables';
    f.Units='normalized';
    f.Position=[0 -0.2 1 1];
    
    
    def_max=max(oo_.SmoothedVariables.def_rate_m);
    def_min=min(oo_.SmoothedVariables.def_rate_m);
    
    subplot(2,2,1);
    plot(Date,oo_.SmoothedVariables.def_rate_m,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Household Default Rate','FontSize',20);
% ylim([def_min def_max]);

    subplot(2,2,2);
    plot(Date,oo_.SmoothedVariables.def_rate_B,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Bank Default Rate','FontSize',20);
% ylim([def_min def_max]);


   subplot(2,2,3);
    plot(Date,oo_.SmoothedVariables.Vs,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Saver Welfare','FontSize',20);

  subplot(2,2,4);
    plot(Date,oo_.SmoothedVariables.Vm,'color','black','lineWidth',5);
 xlim([startDate endDate])
 datetick('x','yy','keeplimits');
title('Borrower Welfare','FontSize',20);


% set(gca,'FontSize',15);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'smoothed_variables','-dpdf');