clear;
clc;
close all;

syms aa bb;
syms xx yy;

eqn(1)= xx+yy- aa;
eqn(2)= xx-2*yy-bb;

result=solve(eqn,[xx yy]);