% clear;clc;close all;
% load welfare_analysis_thr_moments_10e3.mat;
function [welfare_eps_max adhoc_max_eps optimalPara_eps_welfare optimalPara_eps_adhoc welfare_improvement_eps adhoc_improvement_eps]=...
welfare_analysis_eps();

load  welfare_analysis_thr_moments.mat;


%indices for benchmark calibration
Cyphi_bb=1;
phi_Hs_bb=1;
epsilonH1s_bb=1;

%maximization over epsilonH1s only

welfare_eps=reshape(welfare_(1,1,:),[grid_length 1 1]);
welfare_var_eps=reshape(welfare_var_(1,1,:),[grid_length 1 1]);
adhoc_obj_eps=reshape(adhoc_obj_(1,1,:),[grid_length 1 1]);
adhoc_obj_var_eps=reshape(adhoc_obj_var_(1,1,:),[grid_length 1 1]);

objective=weight_level_welfare*welfare_eps-weight_var_welfare*welfare_var_eps;

welfare_eps_max=max(objective);
max_ind=find(objective==welfare_eps_max);
optimalPara_eps_welfare=[epsilonH1s_grid(max_ind)];


objective=weight_level_adhoc*adhoc_obj_eps-weight_var_adhoc*adhoc_obj_var_eps;
adhoc_max_eps=max(objective);
max_ind=find(objective==adhoc_max_eps);
optimalPara_eps_adhoc=[epsilonH1s_grid(max_ind)];

%===============================
welfare_improvement_eps=100*abs((welfare_eps_max-welfare_baseline)/welfare_baseline);

adhoc_improvement_eps=100*(adhoc_max_eps-adhoc_baseline)/adhoc_baseline;
