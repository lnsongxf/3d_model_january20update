clear;clc;close all
%current
ws=1.07;L_ss=0.69;delta_H=0.01;H_ss=12.28;C_ss=1.5;C_bs=0.44;C_es=0.44;a_s=0.5;Trs=0;PIs=0;PHs=0;R_DDs=1.005;
 Ds_current= (ws*L_ss - delta_H*H_ss -C_ss + C_bs + C_es - a_s*Trs+PIs+PHs)/(1-R_DDs)
%original
ws= 1.4167;L_ss= 0.787947;delta_H= 0.01;H_ss=24.7229;C_ss= 1.48586;C_bs=0.0272561;C_es= 0.552446;a_s=0.5;Trs=0.0271452;PIs=0;PHs=0;R_DDs=1.00503;
Ds_original= (ws*L_ss - delta_H*H_ss -C_ss + C_bs + C_es - a_s*Trs+PIs+PHs)/(1-R_DDs)