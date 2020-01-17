clear;clc;%close all
set_parameter_values;
dynare LTV1_noMeasurement noclearall nolog;
welfare_baseline=oo_.mean(134);
welfare_var_baseline=(oo_.var(134,134));

mean_baseline=oo_.mean;
var_baseline=diag(diag(oo_.var));
var_baseline=(var_baseline);%switch to std dev

grid_length=20;





% phi_Hs=0.11;
% epsilonH1s_grid=linspace(0.95,0.7,grid_length);
epsilonH1s_grid=0.86;
% phi_Hs_grid=linspace(0.11,0.15,grid_length);
phi_Hs_grid=linspace(0.05,0.24,grid_length);
% phi_Hs_grid=0.11;
% Cyphi_grid=linspace(0,1,grid_length);

global M_;


var_names={'b_e','b_m','C_m','C_s', 'def_rate_m','def_rate_B','I','Vs','Vm','Y_net_2','welfare_'};
 var_names_plot={'Business Borrowing','Household Borrowing','Consumption of Borrowers','Consumption of Lenders',...
     'Household Default Rate','Bank Default Rate','Investment',...
    'Welfare of Savers','Welfare of Borrowers','Output','Total Welfare'};


 field_names=cellstr(M_.endo_names);
 index=1;
  for jj=1:length(field_names);
 
if true(strcmp(field_names(jj),var_names(index)))
var_indices(index)=jj;
if index<length(var_names)
    index=index+1;
end

end
      
  end
  %======================================


for jj=1:length(phi_Hs_grid)
% for jj=1:length(Cyphi_grid)
    for ee=1:length(epsilonH1s_grid)
% phis=phi_Hs_grid(jj);
phis=0.11;
% phi_Fs=phi_Hs_grid(jj);
phi_Fs=0.11;
phi_Hs=phi_Hs_grid(jj);
% phi_Hs=0.11;

        
%     epsilonH1s=epsilonH1s_grid(ee);%epsilonH1s_grid(ee);
 epsilonH1s=0.86;   
% epsilonF1s=epsilonH1s_grid(ee);
epsilonF1s=0.86;
% Cyphi_H=Cyphi_grid(jj);
% Cyphi_F=Cyphi_grid(jj);
Cyphi_H=0;
Cyphi_F=0;
LTVHrule=0;
LTVFrule=0;
set_parameter_values_policy( Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,phis,epsilonH1s,epsilonF1s,LTVHrule,LTVFrule);
dynare LTV1_noMeasurement.mod noclearall nolog;
SS(jj,ee,:)=oo_.mean;
SS_var(jj,ee,:,:)=oo_.var;
names=M_.endo_names;



    end
end;


sswelfare=SS(:,:,134);
sswelfare_var=SS_var(:,:,134,134);
sswelfare=squeeze(sswelfare);
SS=squeeze(SS);
SS_var=squeeze(SS_var);

SS_var=(SS_var);%switch to standard deviation

 
figure('Name','change in steady state','units','normalized','outerposition',[0 0 1 1]);
for jj=1:length(var_names)
    subplot(4,3,jj);
%     plot(epsilonH1s_grid(1:end),squeeze(SS(:,var_indices(jj))),'color','black');
%     hold on;
%     plot(epsilonH1s_grid(1:end),ones(length(epsilonH1s_grid(1:end)),1)*mean_baseline(var_indices(jj)));
%     xlim([epsilonH1s_grid(end) epsilonH1s_grid(1)]);

    plot(phi_Hs_grid(1:end),squeeze(SS(:,var_indices(jj))),'color','black');
    hold on;
    plot(phi_Hs_grid(1:end),ones(length(phi_Hs_grid(1:end)),1)*mean_baseline(var_indices(jj)));
    xlim([phi_Hs_grid(1) phi_Hs_grid(end)]);
    title(var_names_plot(jj));

%     plot(Cyphi_grid(1:end),squeeze(SS(:,var_indices(jj))),'color','black');
%     hold on;
%     plot(Cyphi_grid(1:end),ones(length(Cyphi_grid(1:end)),1)*mean_baseline(var_indices(jj)));
%     xlim([Cyphi_grid(1) Cyphi_grid(end)]);
%     title(var_names_plot(jj));

end


fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'WA_SCR_housing_level','-dpdf');

%===========================================================
figure('Name','change in variances','units','normalized','outerposition',[0 0 1 1]);
for jj=1:length(var_names)
    subplot(4,3,jj);
%     plot(epsilonH1s_grid(1:end),sqrt(squeeze(SS_var(:,var_indices(jj),var_indices(jj)))),'color','black');
%        hold on;
%       plot(epsilonH1s_grid(1:end),ones(length(epsilonH1s_grid(1:end)),1)*sqrt(var_baseline(var_indices(jj),var_indices(jj))));
%     xlim([epsilonH1s_grid(end) epsilonH1s_grid(1)]);

    plot(phi_Hs_grid(1:end),sqrt(squeeze(SS_var(:,var_indices(jj),var_indices(jj)))),'color','black');
       hold on;
      plot(phi_Hs_grid(1:end),ones(length(phi_Hs_grid(1:end)),1)*sqrt(var_baseline(var_indices(jj),var_indices(jj))));
    xlim([phi_Hs_grid(1) phi_Hs_grid(end)]);    

%     plot(Cyphi_grid(1:end),sqrt(squeeze(SS_var(:,var_indices(jj),var_indices(jj)))),'color','black');
%        hold on;
%       plot(Cyphi_grid(1:end),ones(length(Cyphi_grid(1:end)),1)*sqrt(var_baseline(var_indices(jj),var_indices(jj))));
%     xlim([Cyphi_grid(1) Cyphi_grid(end)]);   
    
    title(var_names_plot(jj));
end


fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'WA_SCR_housing_var','-dpdf');
%==============================welfare improvements with a single policy

 
figure('Name','welfare-level','units','normalized','outerposition',[0 0 1 1]);
 jj=length(var_names)
%     plot(epsilonH1s_grid(1:end),squeeze(SS(:,var_indices(jj))),'color','black');
%     hold on;
%     plot(epsilonH1s_grid(1:end),ones(length(epsilonH1s_grid(1:end)),1)*mean_baseline(var_indices(jj)));
%     xlim([epsilonH1s_grid(end) epsilonH1s_grid(1)]);

    plot(phi_Hs_grid(1:end),squeeze(SS(:,var_indices(jj))),'color','black','lineWidth',5);
    hold on;
    plot(phi_Hs_grid(1:end),ones(length(phi_Hs_grid(1:end)),1)*mean_baseline(var_indices(jj)),'color','red','lineStyle','--','lineWidth',5);
    xlim([phi_Hs_grid(1) phi_Hs_grid(end)]);
    title(['Level of',var_names_plot(jj)]);

%     plot(Cyphi_grid(1:end),squeeze(SS(:,var_indices(jj))),'color','black');
%     hold on;
%     plot(Cyphi_grid(1:end),ones(length(Cyphi_grid(1:end)),1)*mean_baseline(var_indices(jj)));
%     xlim([Cyphi_grid(1) Cyphi_grid(end)]);
%     title(var_names_plot(jj));



fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'welfare_SCR_housing_level','-dpdf');

%===========================================================
figure('Name','welfare-volatility','units','normalized','outerposition',[0 0 1 1]);
jj=length(var_names)
%     plot(epsilonH1s_grid(1:end),sqrt(squeeze(SS_var(:,var_indices(jj),var_indices(jj)))),'color','black');
%        hold on;
%       plot(epsilonH1s_grid(1:end),ones(length(epsilonH1s_grid(1:end)),1)*sqrt(var_baseline(var_indices(jj),var_indices(jj))));
%     xlim([epsilonH1s_grid(end) epsilonH1s_grid(1)]);

    plot(phi_Hs_grid(1:end),sqrt(squeeze(SS_var(:,var_indices(jj),var_indices(jj)))),'color','black','lineWidth',5);
       hold on;
      plot(phi_Hs_grid(1:end),ones(length(phi_Hs_grid(1:end)),1)*sqrt(var_baseline(var_indices(jj),var_indices(jj))),'color','red','lineStyle','--','lineWidth',5);
    xlim([phi_Hs_grid(1) phi_Hs_grid(end)]);    

%     plot(Cyphi_grid(1:end),sqrt(squeeze(SS_var(:,var_indices(jj),var_indices(jj)))),'color','black');
%        hold on;
%       plot(Cyphi_grid(1:end),ones(length(Cyphi_grid(1:end)),1)*sqrt(var_baseline(var_indices(jj),var_indices(jj))));
%     xlim([Cyphi_grid(1) Cyphi_grid(end)]);  

title(['Volatility of ',var_names_plot(jj)]);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'welfare_SCR_housing_var','-dpdf');














var_grid=linspace(0,0.1,5);

for jj=1:length(var_grid)

weight_var_welfare=var_grid(jj);
weight_var_adhoc=0;
weight_level_welfare=1;
weight_level_adhoc=1;

objective_baseline=weight_level_welfare*welfare_baseline-...
    weight_var_welfare*sqrt(welfare_var_baseline);


objective=weight_level_welfare*sswelfare-...
    weight_var_welfare*sqrt(sswelfare_var);

welfare_max=max(objective);
ind_aux=find(objective==welfare_max);
ind_aux=ind_aux(1);
% optimPara_welfare(jj,:)=epsilonH1s_grid(ind_aux)
optimPara_welfare(jj,:)=phi_Hs_grid(ind_aux)
% optimPara_welfare(jj,:)=Cyphi_grid(ind_aux)




welfare_improvement(jj)=100*abs((welfare_max-objective_baseline)/objective_baseline)

end














