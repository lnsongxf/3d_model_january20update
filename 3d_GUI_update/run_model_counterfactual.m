function counterfactual= run_model_counterfactual(data)    
set_parameter_values_policy(data);
dynare LTV1.mod noclearall;
evalin('base','save model_workspace_counterfactual.mat');
load model_workspace_counterfactual.mat;
var_names=who;
for jj=1:length(var_names);
counterfactual.(var_names{jj}) = eval(var_names{jj});
end


end