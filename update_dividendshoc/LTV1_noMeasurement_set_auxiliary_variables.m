function y = LTV1_noMeasurement_set_auxiliary_variables(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

y(135)=y(81)*y(49)*y(15);
y(136)=y(80)*((1-y(113))*y(50)*y(16)-y(52)*(1-params(70))*(1-y(32)));
