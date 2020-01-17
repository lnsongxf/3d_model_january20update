function allvars= retrieve_model    
load model_workspace.mat;
var_names=who;
for jj=1:length(var_names);
allvars.(var_names{jj}) = eval(var_names{jj});
end

end