clear;clc;%close all;
main;
startDate=datenum('01-01-1999');
endDate = datenum('01-12-2016');
T=71;%sample length

Date=linspace(startDate,endDate,T);
%problem: first and second smoother yield different results even without
%changing any parameters. Use the first run of simult_ as the baseline
%version for the time being
%oo_.SmoothedVariables==> smoothed data baseline specification
dy_data_real=oo_.SmoothedVariables.dy_data;

 %variable indices: dy_data=180, dq_H=181,d_b_to_Y=186,
 %int_rate_HH=188,bsp_H=189, b_e=1, b_m=7, C=8,
 
 %var_indices=[1,7,8,12,13,14 171 172];
%  actual_names={'Business Debt','Household Debt','Aggregate Consumption','Def Rate E', 'Def Rate M','Def Rate B','Utility of Borrowing HH','Utility of Lending HH','Average Capital Ratio'...
%      'Output Growth','House Price Growth','Credit to GDP Growth','Mortgage Rate','Household Credit Spread'};
%  actual_names={'Business Debt','Household Debt','Aggregate Consumption','Def Rate E', 'Def Rate M','Def Rate B','Output Growth','House Price Growth','Capital Ratio','Household Interest Rate','Corporate Interest Rate'};
%   actual_names=cellstr(actual_names);

%THIS SHOULD BE IN ORDER OF DECLARATION IN .MOD FILE
%   var_names={'b_e','b_m','C','D','def_rate_e','def_rate_m','def_rate_B','R_D','R_F','R_m','R_tilde_F','R_tilde_H','bsp_F' ,'bsp_H',...           
%    'phib','dy_data','dq_H_data'};

% var_names={'b_e','b_m','C_m','C_s','H_m','H_s','L_m','L_s',                     'dy_data','dq_H_data','dw_data','dc_data'};
% var_names_plot={'Business Borrowing','Household Borrowing','Consumption of Borrowers','Consumption of Lenders',...
% 'Housing owned by Borrowers','housing owned by Lenders','Labor Borrowers','Labor Lenders','Output Growth','House Price Growth','Wage Growth','Consumption Growth'};

% var_names={'b_e','b_m','C_m','C_s',       'R_F','R_H'  ,  'R_tilde_F','R_tilde_H',      'phib',    'dy_data','dq_H_data','dw_data','dc_data'};
% var_names_plot={'Business Borrowing','Household Borrowing','Consumption of Borrowers','Consumption of Lenders',...
% 'Business Lending Rate','Household Lending Rate','Business Portfolio Return','Housing Portfolio Return',...
% 'Bank Capital','Output Growth','House Price Growth','Wage Growth','Consumption Growth'};

%  var_names={'b_e','b_m','C_m','C_s', 'dy_data','dq_H_data','dbe_data','dbm_data'};
%  var_names_plot={'Business Borrowing','Household Borrowing','Consumption of Borrowers','Consumption of Lenders',...
% 'Output Growth','House Price Growth','Business Credit Growth','Mortgage Credit Growth'};

var_names={'b_e','b_m','C_m','C_s', 'def_rate_m','def_rate_B','I','q_H','Vs','Vm','Y_net_2','welfare_'};
 var_names_plot={'Business Borrowing','Household Borrowing','Consumption of Borrowers','Consumption of Lenders',...
     'Household Default Rate','Bank Default Rate','Investment','House Prices',...
    'Welfare of Savers','Welfare of Borrowers','Output','Total Welfare'};

%, ...
%'Output Growth','House Price Growth','Business Credit Growth','Mortgage Credit Growth'

%   var_indices=zeros(length(var_names),1);
% var_names=cellstr(var_names);
%   global M_;
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

%  actual_names=actual_names';
 %extract the shocks from baseline parameterization
ex_BASELINE=[ oo_.SmoothedShocks.epsiA  oo_.SmoothedShocks.epsiJ     oo_.SmoothedShocks.epsiK  ...
    oo_.SmoothedShocks.epsiSe  oo_.SmoothedShocks.epsiSm oo_.SmoothedShocks.epsiSB  ...
    oo_.SmoothedShocks.epsiWb  oo_.SmoothedShocks.epsiWe  ...
    oo_.SmoothedShocks.epsiH oo_.SmoothedShocks.epsiHd  oo_.SmoothedShocks.epsiHk ...
    oo_.SmoothedShocks.epsimarkup_m oo_.SmoothedShocks.epsimarkup_F...
    oo_.SmoothedShocks.epsiEC oo_.SmoothedShocks.epsiECAB oo_.SmoothedShocks.epsiEL ...
    oo_.SmoothedShocks.epsiEbH oo_.SmoothedShocks.epsiEbF...
    oo_.SmoothedShocks.epsi_phiH oo_.SmoothedShocks.epsi_phiF oo_.SmoothedShocks.epsi_epsH];
   
 y0_BASELINE=oo_.dr.ys;
% y0_BASELINE=zeros(length(oo_.dr.ys),1);
dr_BASELINE=oo_.dr;
iorder=1;
%re-simulate the system--> this should be identical/close to the model output
 y_=simult_(y0_BASELINE,dr_BASELINE,ex_BASELINE,iorder);
%y_ ==> variables in the baseline simulation
% variables_BASELINE=[y_(1,1:end);y_(7,1:end);y_(8,1:end);y_(12,1:end);y_(13,1:end);y_(14,1:end);y_(171,1:end);y_(172,1:end)];
variables_BASELINE=y_;
%set prudential policy in the counterfactual exercise
Cyphi_H=0;
Cyphi_F=0;
phi_Hs=0.11;%;0.2112;
phi_Fs=0.11;%0.05;
epsilonH1s=0.86;
epsilonF1s=0.86;
LTVHrule=0;
LTVFrule=0;
phis=0.11;
 zeta_m=0.1;%[0.691226178352873]; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
zeta_F=0.1;%[0.447236412060053]; %SS*** Interest rate stickiness,  appears in 51 FOC for bank business lending, 51 is not in use  %60
 cf_no='11';
% estimated values ---|>%  set_zeta_m=  0.6292
% set_zeta_F=0.7819; 

%add permanent shocks :
%19: phi_H,20: phi_F, 21: epsH--> indices of shocks in ex_BASELINE
%133 phi_F, 134 phi_H, 15 epsilonH --> indices of variables in M_endo_names
%these are not exact, something in log-linearization messes things up

%order: phi_H, phi_F epsilonH
init_period=10;
total=0.0488;duration=20;
%total=0;duration=20;
ex_BASELINE(init_period:init_period+duration-1,19)=(total/duration)*ones(duration,1);%
total=0.015;duration=20;
%total=0;duration=20;
ex_BASELINE(init_period:init_period+duration-1,20)=(total/duration)*ones(duration,1);
total=0.0808;duration=20;
%total=0;duration=20;
ex_BASELINE(init_period:init_period+duration-1,21)=(total/duration)*ones(duration,1);

%cf9=zeta_m=0.1,phi_Hs=0.15; cf10=corporate stickiness
%=================================================
% set_parameter_values_zeta(zeta_m,zeta_F);
%set_parameter_values_policy(Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,phis,epsilonH1s,epsilonF1s,LTVHrule,LTVFrule);
set_parameter_values_counterfactual(Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,phis,epsilonH1s,epsilonF1s,LTVHrule,LTVFrule,zeta_m,zeta_F);
dynare LTV1.mod noclearall nolog;
y0_CF=oo_.dr.ys;
dr_CF=oo_.dr;
iorder=1;
 y_CF=simult_(y0_BASELINE,dr_CF,ex_BASELINE,iorder);%smoother swith baseline initial values and shocks. decision rules are different
%y_ ==> smoothed data counterfactual specification
% variables_CF=[y_CF(1,1:end);y_CF(7,1:end);y_CF(8,1:end);y_CF(12,1:end);y_CF(13,1:end);y_CF(14,1:end);y_CF(171,1:end);y_CF(172,1:end)];
variables_CF=y_CF;
    f=figure;
    f.Name='Counterfactuals';
    f.Units='normalized';
    f.Position=[0 -0.2 1 1];
    index=0;
    
for jj=var_indices
    index=index+1;
subplot(4,3,index);

plot(Date,variables_BASELINE(jj,3:end),'color','black','lineWidth',1);
hold on;
plot(Date,variables_CF(jj,3:end),'--','color','red','lineWidth',1);

% if jj==var_indices(end)
if jj==1
legend('baseline','counterfactual');
end
  xlim([startDate endDate])
  datetick('x','yy','keeplimits');
%       title(field_names(jj));
title(var_names_plot(index),'FontSize',15);
%       set(gca,'FontSize',13);
end
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
string=['counterfactuals' cf_no];
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,string,'-dpdf');


%summary statistics
round_dec=4;
%baseline_mean=round(mean(variables_BASELINE(:,3:end)')',4);
%baseline_var=round(sqrt(var(variables_BASELINE(:,3:end)'))',4);%switch to std. dev
%cf_mean=round(mean(variables_CF(:,3:end)')',4);
%cf_var=round(sqrt(var(variables_CF(:,3:end)'))',4);%switch to std. dev
baseline_mean=mean(variables_BASELINE(:,3:end)')';
baseline_var=sqrt(var(variables_BASELINE(:,3:end)'))';%switch to std. dev
cf_mean=mean(variables_CF(:,3:end)')';
cf_var=sqrt(var(variables_CF(:,3:end)'))';%switch to std. dev
perc_change_mean=(cf_mean-baseline_mean)./baseline_mean;
perc_change_var=(cf_var-baseline_var)./baseline_var;
%=========================
field_names=field_names(var_indices);
baseline_mean=baseline_mean(var_indices);
baseline_var=baseline_var(var_indices);
cf_mean=cf_mean(var_indices);
cf_var=cf_var(var_indices);
perc_change_mean=perc_change_mean(var_indices);
perc_change_var=perc_change_var(var_indices);

summary_stat=table(field_names,baseline_mean,cf_mean, perc_change_mean,baseline_var,cf_var, perc_change_var)


    f=figure;
    f.Name='Counterfactuals II';
    f.Units='normalized';
    f.Position=[0 0.2 1 0.5];
    index=0;
    

subplot(1,3,1);
plot(Date,variables_BASELINE(134,3:end),'color','black','lineWidth',1);
hold on;
plot(Date,variables_CF(134,3:end),'--','color','red','lineWidth',1);
title('Sectoral Capital Requirements: Housing');
  xlim([startDate endDate])
  datetick('x','yy','keeplimits');

subplot(1,3,2);
plot(Date,variables_BASELINE(133,3:end),'color','black','lineWidth',1);
hold on;
plot(Date,variables_CF(133,3:end),'--','color','red','lineWidth',1);
title('Sectoral Capital Requirements: Corporate');
  xlim([startDate endDate])
  datetick('x','yy','keeplimits');

subplot(1,3,3);
plot(Date,variables_BASELINE(15,3:end),'color','black','lineWidth',1);
hold on;
plot(Date,variables_CF(15,3:end),'--','color','red','lineWidth',1);
title('LTV Limit');
  xlim([startDate endDate])
  datetick('x','yy','keeplimits');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
string=['CF_policy_rules' cf_no];
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,string,'-dpdf');


% figure;
% plot(dy_data_real(2:end));
% hold on;
% plot(y_CF(180,3:end));

% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'counterfactual_lowLTV','-dpdf');