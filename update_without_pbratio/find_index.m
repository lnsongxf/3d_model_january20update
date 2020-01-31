function [ind]=find_index(number,grid_length);

% ind1=ceil(number/grid_length^2);
% 
% number2=   (number-(ind1-1)*grid_length^2);
% ind2=ceil(number2/grid_length);
% 
% number3= number2-(ind2-1)*grid_length;
% ind3=number3;
% 
% ind=[ind1,ind2,ind3];



ind3=ceil(number/grid_length^2);

number2=   (number-(ind3-1)*grid_length^2);
ind2=ceil(number2/grid_length);

number3= number2-(ind2-1)*grid_length;
ind1=number3;

ind=[ind1,ind2,ind3];

end