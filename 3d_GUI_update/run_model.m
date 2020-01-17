function allvars= run_model(data)    
set_parameter_values(data);
dynare LTV1.mod noclearall;



evalin('base','save model_workspace.mat '); % -regexp ^(?!(ff|ans)$).
load model_workspace.mat;
clearvars ff;

var_names=who;
for jj=1:length(var_names);
allvars.(var_names{jj}) = eval(var_names{jj});
end

% allvars.M_= eval('M_');
% allvars.options_= eval('options_');
%  why do these not end up in the first list? 


end