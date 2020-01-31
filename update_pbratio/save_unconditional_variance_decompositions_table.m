var_decomp=round(oo_.variance_decomposition,2); %size of numEndo x numShock


%var_names={'b_e','b_m','C_m','C_s', 'def_rate_m','def_rate_B','I','q_H','R_m','Vs','Vm','Y_net_2'};
 %var_names_plot={'Business Borrowing','Household Borrowing','Consumption of Borrowers','Consumption of Lenders',...
 %    'Household Default Rate','Bank Default Rate','Investment','House Prices','Mortgage Rate',...
  %  'Welfare of Savers','Welfare of Borrowers','Output'};
  
var_names={'def_rate_m','def_rate_B','R_D','R_F','R_m','Vs','Vm','dy_data','dq_H_data','dbe_data','dbm_data','dw_data','dinve_data','dc_data'};
  
shock_names={'A','J','K','Se','Sm','SB','Wb','We','H','Hd','Hk','Markup m','Markup e','EC','ECAB'};
shock_names=shock_names';
%, ...
%'Output Growth','House Price Growth','Business Credit Growth','Mortgage Credit Growth'

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
  
  
  output_file='variance_decompositions.csv';
 output_sheet='variance_decompositions';
  xlswrite(output_file,var_names,output_sheet,'b1');
  xlswrite(output_file,shock_names,output_sheet,'a2');
   xlswrite(output_file,var_decomp(var_indices,:)',output_sheet,'b2');
  
  