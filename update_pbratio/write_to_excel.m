estim_para_names2=(fieldnames(oo_.posterior_mode.parameters));
estim_para_names1=(fieldnames(oo_.posterior_mode.shocks_std));
estim_para_names=[estim_para_names1; estim_para_names2];

estim_para_mode=oo_.posterior.optimization.mode;
estim_para_std=sqrt(diag(oo_.posterior.optimization.Variance));

% laplace=NaN;
% laplace=oo_.posterior.optimization.log_density;

  output_file='Estimation_Results.csv';
 output_sheet='Estimation_Results';
 
  estim_para_mode=round(estim_para_mode,2);
estim_para_std=round(estim_para_std,3);

 xlswrite(output_file,estim_para_names,output_sheet,'a2');
  xlswrite(output_file,estim_para_mode,output_sheet,'b2');
   xlswrite(output_file,estim_para_std,output_sheet,'c2');
%    
%      xlswrite(output_file,laplace,output_sheet,'d2');