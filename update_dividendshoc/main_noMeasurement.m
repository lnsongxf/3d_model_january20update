clear;clc;%close all;
%==========================================
% SET PARAMETER VALUES========================
%==========================================
set_parameter_values;
% LTV1_steadystate;
dynare LTV1_noMeasurement.mod noclearall;
%=====================================




% variable_names=M_.endo_names;% e.g: b_e=1, b_m=7,C=8
% steady_state=oo_.steady_state;
% 
% table([variable_names(1),variable_names(7),variable_names(8)]',...
% [steady_state(1),steady_state(7),steady_state(8)]')
% steady;