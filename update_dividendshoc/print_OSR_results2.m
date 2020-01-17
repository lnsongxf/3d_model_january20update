clear;clc;close all;
load optimal_simple_rules_housing;
var_grid=linspace(0,0.1,5);
optimalPara_welfare=zeros(3,length(var_grid));
welfare_improvement=zeros(length(var_grid),1);
welfare_baseline=oo_.mean(134);
welfare_var_baseline=sqrt(oo_.var(134,134));
for jj=1:length(var_grid)

weight_var_welfare=var_grid(jj);
weight_var_adhoc=0;
weight_level_welfare=1;
weight_level_adhoc=1;

welfare_baseline1=weight_level_welfare*welfare_baseline-...
    weight_var_welfare*welfare_var_baseline;

% welfare_=(variables_mean(:,:,:,9).*variables_mean(:,:,:,84)+variables_mean(:,:,:,10).*variables_mean(:,:,:,83))./...
%    ( variables_mean(:,:,:,9)+variables_mean(:,:,:,10));
welfare_=variables_mean(:,:,:,134);

% welfare_var_=(variables_mean(:,:,:,9).*variables_var(:,:,:,84)+variables_var(:,:,:,10).*variables_var(:,:,:,83))./...
%    ( variables_var(:,:,:,9)+variables_var(:,:,:,10));

welfare_var_=sqrt(variables_var(:,:,:,134,134));


% adhoc_obj_=variables_mean(:,:,:,1)*weights_adhoc(1)+...
%     variables_mean(:,:,:,7)*weights_adhoc(2)+...
%     variables_mean(:,:,:,109)*weights_adhoc(3);

% adhoc_obj_var_=variables_var(:,:,:,1)*weights_adhoc(1)+...
%     variables_var(:,:,:,7)*weights_adhoc(2)+...
%     variables_var(:,:,:,109)*weights_adhoc(3);
%=======================================================
objective=weight_level_welfare*welfare_-...
    weight_var_welfare*welfare_var_;

welfare_max=max(max(max(objective)));
ind_aux=find(objective==welfare_max);
max_ind=find_index(ind_aux,grid_length);

optimalPara_welfare(:,jj)=[phi_Hs_grid(max_ind(1))...
    ,Cyphi_H_grid(max_ind(2)),epsilonH1s_grid(max_ind(3))];

% objective=weight_level_adhoc*adhoc_obj_-weight_var_adhoc*adhoc_obj_var_;
% adhoc_max=max(max(max(objective)));
% ind_aux=find(objective==adhoc_max);
% max_ind=find_index(ind_aux,grid_length);
% 
% optimalPara_adhoc=[phi_Hs_grid(max_ind(1))...
%     ,phi_Fs_grid(max_ind(2)),epsilonH1s_grid(max_ind(3))];

%===============================
welfare_improvement(jj)=100*abs((welfare_max-welfare_baseline1)/welfare_baseline1);

% adhoc_improvement=100*(adhoc_max-adhoc_baseline)/adhoc_baseline;




end
% disp('AD HOC RULE MAXIMIZING PARAMETERS:')
% disp(optimalPara_adhoc);
% 
% disp('PERCENTAGE IMPROVEMENT IN WELFARE COMPARED TO BASELINE:')
% disp(adhoc_improvement);
disp('WELFARE MAXIMIZING PARAMETERS:')
disp(optimalPara_welfare);

disp('PERCENTAGE IMPROVEMENT IN WELFARE COMPARED TO BASELINE:')
disp(welfare_improvement);
% save optimal_simple_rules_housing;


%  [welfare_eps_max adhoc_max_eps optimalPara_eps_welfare optimalPara_eps_adhoc welfare_improvement_eps adhoc_improvement_eps]=...
%      welfare_analysis_eps();
% 
%   [welfare_phi_max adhoc_max_phi optimalPara_phi_welfare optimalPara_phi_adhoc welfare_improvement_phi adhoc_improvement_phi]=...
%      welfare_analysis_phi();
