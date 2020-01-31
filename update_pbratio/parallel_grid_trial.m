clear;clc;close all;


xx=nan(10,10,10);


for jj=1:10
    for ii=1:10
        parfor kk=1:10
           dynare LTV1.mod noclearall nolog;
           
        end
    end
end