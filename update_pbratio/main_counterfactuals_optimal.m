clear;clc;%close all;
main;
startDate=datenum('01-01-1998');
endDate = datenum('01-12-2016');
T=71;%sample length
iorder=1;

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
    oo_.SmoothedShocks.epsiEbH oo_.SmoothedShocks.epsiEbF];
   
 y0_BASELINE=oo_.dr.ys;
% y0_BASELINE=zeros(length(oo_.dr.ys),1);
dr_BASELINE=oo_.dr;
iorder=1;
%re-simulate the system--> this should be identical/close to the model output
 y_=simult_(y0_BASELINE,dr_BASELINE,ex_BASELINE,iorder);
%y_ ==> variables in the baseline simulation
% variables_BASELINE=[y_(1,1:end);y_(7,1:end);y_(8,1:end);y_(12,1:end);y_(13,1:end);y_(14,1:end);y_(171,1:end);y_(172,1:end)];
variables_BASELINE=y_;
baseline_mean=round(mean(variables_BASELINE(:,3:end)')',4);
baseline_var=round(var(variables_BASELINE(:,3:end)')',4);
baseline_mean=baseline_mean(var_indices);
baseline_var=baseline_var(var_indices);

welfare_level_baseline=baseline_mean(end);
welfare_var_baseline=baseline_var(end);
%compute baseline paths up to this point
%=============================================================================================\
%start of counterfactuals
%set prudential policy in the counterfactual exercise

grid_length=10;
phi_Hs_grid=linspace(0.05,0.25,grid_length);
Cyphi_grid=linspace(0,0.5,grid_length);
epsilon_H1s_grid=linspace(0.95,0.7,grid_length);


for jj=1:length(Cyphi_grid)
%    for ee=1:length(epsilon_H1s_grid)
%        for kk=1:length(phi_Hs_grid)
    
% Cyphi_H=0;
Cyphi_H=Cyphi_grid(jj);
% Cyphi_F=0;
% Cyphi_F=Cyphi_grid(jj);
phi_Hs=0.11;
phi_Fs=0.11;
epsilonH1s=0.86*1*1;
epsilonF1s=0.86*1*1;
LTVHrule=0;
LTVFrule=0;
phis=0.11;
 zeta_m=0.05;%[0.691226178352873]; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
zeta_F=0.05;%[0.447236412060053];
% estimated values ---|>%  set_zeta_m=  0.6292
% set_zeta_F=0.7819; 
cf_no='1';
%=================================================
% set_parameter_values_zeta(zeta_m,zeta_F);
%set_parameter_values_policy(Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,phis,epsilonH1s,epsilonF1s,LTVHrule,LTVFrule);
set_parameter_values_counterfactual(Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,phis,epsilonH1s,epsilonF1s,LTVHrule,LTVFrule,zeta_m,zeta_F);
dynare LTV1.mod noclearall nolog;
y0_CF=oo_.dr.ys;
dr_CF=oo_.dr;

 y_CF(:,:,jj)=simult_(y0_BASELINE,dr_CF,ex_BASELINE,iorder);
 variables_CF=y_CF(:,:,jj);
cf_mean=round(mean(variables_CF(:,3:end)')',4);
cf_var=round(var(variables_CF(:,3:end)')',4);
%=========================


cf_mean_grid(:,jj)=cf_mean(var_indices);
cf_var_grid(:,jj)=cf_var(var_indices);


welfare_level(jj)=cf_mean_grid(end,jj);
welfare_var(jj)=cf_var_grid(end,jj);

end




var_grid=linspace(0,1,10);

for jj=1:length(var_grid)

weight_var_welfare=var_grid(jj);
weight_var_adhoc=0;
weight_level_welfare=1;
weight_level_adhoc=1;

objective_baseline=weight_level_welfare*welfare_level_baseline-...
    weight_var_welfare*sqrt(welfare_var_baseline);


objective=weight_level_welfare*welfare_level-...
    weight_var_welfare*sqrt(welfare_var);

welfare_max=max(objective);
ind_aux=find(objective==welfare_max);
ind_aux=ind_aux(1);
% optimPara_welfare(jj,:)=epsilonH1s_grid(ind_aux)
% optimPara_welfare(jj,:)=phi_Hs_grid(ind_aux)
optimPara_welfare(jj,:)=Cyphi_grid(ind_aux)




welfare_improvement(jj)=100*abs((welfare_max-objective_baseline)/objective_baseline)

end


