clear;clc;close all;

%policy parameters: phi_Hs,phi_Fs,Cyphi_Hs,Cyphi_Fs, epsilon_H1s
%target variables of interest:
%1) ad-hoc rule: b_e, b_m, Y_net--> indices: 1, 7, 136
%2) welfare-based rule: Vm 110, Vs 109, Cm 9, Cs 10 

%1-shot simulation stored in: oo_.endo_simul
%steady state values stored in oo_.steady_state
%smoothed parameters stored in oo_.SmoothedVariables
%theoretical mean values stored in: oo_.mean
%theoretical variances stored in: oo_.variance
%--> these are returning NaN at the moment

grid_length=5;
weights_adhoc=[1/3,1/3,1/3];
numVar=188;%length(oo_.steady_state)

Cyphi_H_grid=linspace(0,0,grid_length);
Cyphi_F_grid=linspace(0,0,grid_length);
% Cyphi_H_grid=0;
% Cyphi_F_grid=0;
phi_Hs_grid=linspace(0.11,0.2,grid_length);
phi_Fs_grid=linspace(0.11,0.2,grid_length);
epsilonH1s_grid=0.86;

%first try: no sectoral requirements
for cc_ind=1:length(Cyphi_H_grid)

%     Cyphi_H=0;
%     Cyphi_F=0;
         Cyphi_H=Cyphi_H_grid(cc_ind);
         Cyphi_F=Cyphi_F_grid(cc_ind);
    
    for pp_ind=1:length(phi_Hs_grid)
        
        string=['grid no: ' num2str((cc_ind-1)*grid_length+pp_ind)]
        disp(string)
        
        

           phi_Hs=phi_Hs_grid(pp_ind);
           phi_Fs=phi_Fs_grid(pp_ind);
           epsilonH1s=epsilonH1s_grid;
           
set_parameter_values_OSR(Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,epsilonH1s);
try
dynare LTV1.mod noclearall;
simul_mean(:,cc_ind,pp_ind)=mean(oo_.endo_simul');
simul_var(:,cc_ind,pp_ind)=var(oo_.endo_simul');
catch
simul_mean(:,cc_ind,pp_ind)=-Inf*ones(numVar,1);
simul_var(:,cc_ind,pp_ind)=-Inf*ones(numVar,1);
end

% steady_state(steady_state==NaN)=-Inf;
% simul_mean(simul_mean==NaN)=-Inf;
% simul_var(simul_var==NaN)=-Inf;


    end
end

%1) ad-hoc rule: b_e, b_m, Y_net--> indices: 1, 7, 136
%2) welfare-based rule: Vm 110, Vs 109, Cm 9, Cs 10 

matrix_b_e=simul_mean(1,:,:);
matrix_b_e=squeeze(matrix_b_e);

matrix_b_m=simul_mean(7,:,:);
matrix_b_m=squeeze(matrix_b_m);

matrix_Y_net=simul_mean(136,:,:);
matrix_Y_net=squeeze(matrix_Y_net);

matrix_Vs=simul_mean(109,:,:);
matrix_Vs=squeeze(matrix_Vs);
matrix_Vm=simul_mean(110,:,:);
matrix_Vm=squeeze(matrix_Vm);

welfare_diff=abs(simul_mean(110,:,:)-simul_mean(109,:,:));

%welfare_ss=(steady_state(110,:,:).*steady_state(9,:,:)+steady_state(109,:,:).*steady_state(10,:,:))./...
  %  (steady_state(9,:,:)+steady_state(10,:,:));

welfare_simul_mean=(simul_mean(110,:,:).*simul_mean(9,:,:)+simul_mean(109,:,:).*simul_mean(10,:,:))./...
    (simul_mean(9,:,:)+simul_mean(10,:,:));

welfare_simul_var=(simul_var(110,:,:).*simul_var(9,:,:)+simul_var(109,:,:).*simul_var(10,:,:))./...
    (simul_var(9,:,:)+simul_var(10,:,:));

%adhoc_objective_ss=steady_state(1,:,:)*weights_adhoc(1)+...
  %  steady_state(7,:,:)*weights_adhoc(2)+...
    %steady_state(136,:,:)*weights_adhoc(3);

adhoc_objective_simul_mean=simul_mean(1,:,:)*weights_adhoc(1)+...
    simul_mean(7,:,:)*weights_adhoc(2)+...
    simul_mean(136,:,:)*weights_adhoc(3);

adhoc_objective_simul_var=simul_var(1,:,:)*weights_adhoc(1)+...
    simul_var(7,:,:)*weights_adhoc(2)+...
    simul_var(136,:,:)*weights_adhoc(3);



%welfare_ss=squeeze(welfare_ss);
welfare_simul_mean=squeeze(welfare_simul_mean);
welfare_simul_var=squeeze(welfare_simul_var);
%adhoc_objective_ss=squeeze(adhoc_objective_ss);
adhoc_objective_simul_mean=squeeze(adhoc_objective_simul_mean);
adhoc_objective_simul_var=squeeze(adhoc_objective_simul_var);
%welfare_diff=squeeze(welfare_diff);
%welfare2_simul_mean=squeeze(welfare_simul_mean-10e10*welfare_diff);

%welfare_ss_max=max(max(welfare_ss));
%[max_ind(1),max_ind(2)]=find(welfare_ss==welfare_ss_max);

%optimalPara_welfare_ss=[Cyphi_H_grid(max_ind(1))...
  %  ,phi_Hs_grid(max_ind(2))];
%=======================================
welfare_simul_max=max(max(welfare_simul_mean));
[max_ind(1),max_ind(2)]=find(welfare_simul_mean==welfare_simul_max);

optimalPara_welfare_simul=[Cyphi_H_grid(max_ind(1))...
    ,phi_Hs_grid(max_ind(2))];
%============================================
%welfare2_simul_max=max(max(welfare2_simul_mean));
%[max_ind(1),max_ind(2)]=find(welfare2_simul_mean==welfare2_simul_max);

%optimalPara_welfare2_simul=[Cyphi_H_grid(max_ind(1))...
 %   ,phi_Hs_grid(max_ind(2))];

%=======================================
%adhoc_objective_ss_max=max(max(adhoc_objective_ss));
%[max_ind(1),max_ind(2)]=find(adhoc_objective_ss==adhoc_objective_ss_max);

%optimalPara_adhoc_ss=[Cyphi_H_grid(max_ind(1))...
 %   ,phi_Hs_grid(max_ind(2))];
%=======================================
adhoc_objective_simul_max=max(max(adhoc_objective_simul_mean));
[max_ind(1),max_ind(2)]=find(adhoc_objective_simul_mean==adhoc_objective_simul_max);

optimalPara_adhoc_mean=[Cyphi_H_grid(max_ind(1))...
    ,phi_Hs_grid(max_ind(2))];
%=======================================
%welfare_savers_max=max(max(matrix_Vs));
%[max_ind(1),max_ind(2)]=find(matrix_Vs==welfare_savers_max);

%optimalPara_welfare_savers=[Cyphi_H_grid(max_ind(1))...
 %   ,phi_Hs_grid(max_ind(2))];

%welfare_borrowers_max=max(max(matrix_Vm));
%[max_ind(1),max_ind(2)]=find(matrix_Vm==welfare_borrowers_max);

%optimalPara_welfare_borrowers=[Cyphi_H_grid(max_ind(1))...
  %  ,phi_Hs_grid(max_ind(2))];



disp('WELFARE MAXIMIZING PARAMETERS FROM SIMULATION:')
disp(optimalPara_welfare_simul);

disp('AD HOC RULE MAXIMIZING PARAMETERS FROM SIMULATION:')
disp(optimalPara_adhoc_mean);

%disp('WELFARE SAVERS FROM SIMULATION:')
%disp(optimalPara_welfare_savers);

%disp('WELFARE BORROWERS PARAMETERS FROM SIMULATION:')
%disp(optimalPara_welfare_borrowers);


%disp('WELFARE RULE (INCLUDING INEQUALITY) MAXIMIZING PARAMETERS FROM SIMULATION:')
%disp(optimalPara_welfare2_simul);






