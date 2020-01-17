clear;clc;
rng(1);
tic
numPeriods=10000;
drop=5000;
numSimul=100;
numShocks=18;
numVar=188;
iorder=1;%perturbation order in simulations

sigma1=[0.002342077555948;%A
    0.574132008875862;%J
    0;%K
    0.049174575223881;%Se
    0.048120855178746;%Sm
    0.026676974833692%SB
    0;%Wb;
    0.006883356991601;%;We
    0;%H
    0.009266823021518;%Hd
    0.003222749799409;%Hk
    4.491303201735505e-04;%epsimarkup_m
    4.565761897299154e-04;%epsimarkup_F
    0.015487379825398;%EC
     0.069847979400518;%ECAB
     0;%EL
     0;%EbH
     0]; %EbF

var_names={'b_e','b_m','C_m','C_s', 'def_rate_m','def_rate_B','Vs', 'Vm' 'Y_net', 'phib'};
 var_names_plot={'Business Borrowing','Household Borrowing','Consumption of Borrowers','Consumption of Lenders',...
'Household Default Rate','Bank Default Rate',...
 'Welfare of Borrowers','Welfare of Savers','Output','Bank Capital'};


  global M_;
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

exo_names=M_.exo_names;
%exo variables order: 
%epsiA       % epsiJ       % epsiK       % epsiSe      % epsiSm      
% epsiSB      % epsiWb      % epsiWe      % epsiH       % epsiHd      % epsiHk      % epsimarkup_m% epsimarkup_F
% epsiEC      % epsiECAB    % epsiEL      % epsiEbH     % epsiEbF     

%std error vector

   %fixed shocks across all simulations  
%    simul_shocks=zeros(numPeriods,numShocks,numSimul);
% for jj=1:numSimul
%      disp('simulating shocks')
%     disp(jj)
% simul_shocks(:,:,jj)=mvnrnd(zeros(length(sigma1),1),diag(sigma1),numPeriods);
% end

% save simul_shocks.mat simul_shocks;
load simul_shocks.mat simul_shocks;


%=============retrieve the decision rules oo_.dr
Cyphi_H=0;
Cyphi_F=0;
phi_Hs=0.11;
phi_Fs=0.11;
epsilonH1s=0.86;
%============
set_parameter_values_OSR(Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,epsilonH1s);
dynare LTV1.mod noclearall;
dr_current=oo_.dr;
% y0_current=oo_.dr.ys;
y0_current=zeros(numVar,1);

%=======================================


%   simul_shocks==>zeros(numPeriods,numShocks,numSimul);
y_current=zeros(numVar,numPeriods+1,numSimul);
for jj=1:numSimul
    shock_current=simul_shocks(:,:,jj)';
    disp(jj)
    y_current(:,:,jj)=simult_(y0_current,dr_current,simul_shocks(:,:,jj),iorder);
    
end
    


toc
