%use no-measurement file--> variable indices will differ
%no sectoral requirements here yet--

clear;clc;%close all
tic
grid_length=5;
weights_adhoc=[1/3,1/3,1/3];

set_parameter_values;
dynare LTV1_noMeasurement noclearall nolog;

welfare_baseline=(oo_.mean(9)*oo_.mean(84)+oo_.mean(10)*oo_.mean(83))./...
    (oo_.mean(9)+oo_.mean(10));


adhoc_baseline=oo_.mean(1)*weights_adhoc(1)+...
    oo_.mean(7)*weights_adhoc(2)+...
    oo_.mean(109)*weights_adhoc(3);


% phi_Hs=0.11;
% epsilonH1s_grid=0.86;
% phi_Hs_grid=linspace(0.11,0.15,grid_length);
% phi_Hs_grid=linspace(0.11,0.19,grid_length);

global M_;

% var_names={'b_e','b_m','C_m','C_s',  'H_m','H_s','L_m','L_s', 'R_F','R_H'  , 'Vs', 'Vm' 'Y_net', 'phib'};
% var_names_plot={'Business Borrowing','Household Borrowing','Consumption of Borrowers','Consumption of Lenders',...
% 'Housing Owned by Borrowers','Housing Owned by Lenders','Borrowers Labor','Lenders Labor','Business Lending Rate','Household Lending Rate'...
% 'Welfare of Borrowers','Welfare of Savers','Output','Bank Capital'};

 var_names={'b_e','b_m','C_m','C_s', 'def_rate_m','def_rate_B','Vs', 'Vm' 'Y_net_2', 'phib'};
 var_names_plot={'Business Borrowing','Household Borrowing','Consumption of Borrowers','Consumption of Lenders',...
'Household Default Rate','Bank Default Rate',...
 'Welfare of Borrowers','Welfare of Savers','Output','Bank Capital'};

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




numVar=length(oo_.mean);

% Cyphi_H_grid=linspace(0,0,grid_length);
% Cyphi_F_grid=linspace(0,0,grid_length);
Cyphi_H_grid=0;
Cyphi_F_grid=0;
% Cyphi_H_grid=0;
% Cyphi_F_grid=0;
phi_Hs_grid=linspace(0.11,0.2,grid_length);
phi_Fs_grid=linspace(0.11,0.2,grid_length);
epsilonH1s_grid=linspace(0.86,0.7,grid_length);
% epsilon_H1s_grid=[0.86,0.76];
% epsilonH1s_grid=0.86;

%indices of benchmark values




%first try: no sectoral requirements
index=0;
for cc_ind=1:length(Cyphi_H_grid)
for pp_ind=1:length(phi_Hs_grid)
    for ee_ind=1:length(epsilonH1s_grid)
        index=index+1
        toc
        
        string=['grid no: ' num2str(index)]
        disp(string)
        
        
         Cyphi_H=Cyphi_H_grid(cc_ind);
         Cyphi_F=Cyphi_F_grid(cc_ind);
           phi_Hs=phi_Hs_grid(pp_ind);
           phi_Fs=phi_Fs_grid(pp_ind);
           epsilonH1s=epsilonH1s_grid(ee_ind);
           
set_parameter_values_OSR(Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,epsilonH1s);
try
dynare LTV1_noMeasurement.mod noclearall nolog ;
variables_mean(cc_ind,pp_ind,ee_ind,:)=oo_.mean;
variables_var(cc_ind,pp_ind,ee_ind,:,:)=oo_.var;
catch
variables_mean(cc_ind,pp_ind,ee_ind,:)=-Inf*ones(numVar,1);
variables_var(cc_ind,pp_ind,ee_ind,:,:)=-Inf*ones(numVar,1);
end




    end
end
end

%display target variable means
% figure;
% for jj=1:length(var_names);
%     subplot(4,3,jj);
%     surf(variables_mean(:,:,var_indices(jj)));
%     title(var_names(jj))
% end

%====================welfare measure
%indices: C_m 9, C_s 10, V_s 83, V_m 84 
%b_e 1 b_m 7 Y_net 109


weight_var_welfare=0;
weight_var_adhoc=0;
weight_level_welfare=1;
weight_level_adhoc=1;



welfare_=(variables_mean(:,:,:,9).*variables_mean(:,:,:,84)+variables_mean(:,:,:,10).*variables_mean(:,:,:,83))./...
   ( variables_mean(:,:,:,9)+variables_mean(:,:,:,10));

welfare_var_=(variables_var(:,:,:,9).*variables_var(:,:,:,84)+variables_var(:,:,:,10).*variables_var(:,:,:,83))./...
   ( variables_var(:,:,:,9)+variables_var(:,:,:,10));

adhoc_obj_=variables_mean(:,:,:,1)*weights_adhoc(1)+...
    variables_mean(:,:,:,7)*weights_adhoc(2)+...
    variables_mean(:,:,:,109)*weights_adhoc(3);

adhoc_obj_var_=variables_var(:,:,:,1)*weights_adhoc(1)+...
    variables_var(:,:,:,7)*weights_adhoc(2)+...
    variables_var(:,:,:,109)*weights_adhoc(3);
%=======================================================
objective=weight_level_welfare*welfare_-weight_var_welfare*welfare_var_;

welfare_max=max(max(max(objective)));
ind_aux=find(objective==welfare_max);
max_ind=find_index(ind_aux,grid_length);

optimalPara_welfare=[Cyphi_H_grid(max_ind(1))...
    ,phi_Hs_grid(max_ind(2)),epsilonH1s_grid(max_ind(3))];

objective=weight_level_adhoc*adhoc_obj_-weight_var_adhoc*adhoc_obj_var_;
adhoc_max=max(max(max(objective)));
ind_aux=find(objective==adhoc_max);
max_ind=find_index(ind_aux,grid_length);

optimalPara_adhoc=[Cyphi_H_grid(max_ind(1))...
    ,phi_Hs_grid(max_ind(2)),epsilonH1s_grid(max_ind(3))];

%===============================
welfare_improvement=100*abs((welfare_max-welfare_baseline)/welfare_baseline);

adhoc_improvement=100*(adhoc_max-adhoc_baseline)/adhoc_baseline;


disp('WELFARE MAXIMIZING PARAMETERS:')
disp(optimalPara_welfare);

disp('PERCENTAGE IMPROVEMENT IN WELFARE COMPARED TO BASELINE:')
disp(welfare_improvement);


disp('AD HOC RULE MAXIMIZING PARAMETERS:')
disp(optimalPara_adhoc);

disp('PERCENTAGE IMPROVEMENT IN WELFARE COMPARED TO BASELINE:')
disp(adhoc_improvement);

save welfare_analysis_thr_moments.mat;


 [welfare_eps_max adhoc_max_eps optimalPara_eps_welfare optimalPara_eps_adhoc welfare_improvement_eps adhoc_improvement_eps]=...
     welfare_analysis_eps();

  [welfare_phi_max adhoc_max_phi optimalPara_phi_welfare optimalPara_phi_adhoc welfare_improvement_phi adhoc_improvement_phi]=...
     welfare_analysis_phi();
 

 
 









%optimal rules and improvements--1 at a time







