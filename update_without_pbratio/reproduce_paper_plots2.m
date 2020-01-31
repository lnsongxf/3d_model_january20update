%TURN OFF zeta_m and zeta_F in estimation block of .mod file before running

clear;clc;%close all;
main;

%shock order : epsiA epsiJ epsiSe epsiSm epsiSB epsiWe epsiWb epsiHd epsiHk 
%variable order: b_m b_e R_m R_F Y_net C

%irf order: irf_length, shock, variable
irf_length=40;
numVar=7;
numShock=9;
irf=nan(irf_length,numShock,numVar);

% 
% irf1(:,:,1)=[b_m_epsiA b_m_epsiJ b_m_epsiSm   b_m_epsiWe  b_m_epsiHd b_m_epsiHk b_m_epsiEC b_m_epsiECAB] ;
% irf1(:,:,2)=[b_e_epsiA b_e_epsiJ b_e_epsiSm   b_e_epsiWe  b_e_epsiHd b_e_epsiHk b_e_epsiEC b_e_epsiECAB] ;
% irf1(:,:,3)=[R_m_epsiA R_m_epsiJ R_m_epsiSm   R_m_epsiWe  R_m_epsiHd R_m_epsiHk R_m_epsiEC R_m_epsiECAB] ;
% irf1(:,:,4)=[R_F_epsiA R_F_epsiJ R_F_epsiSm   R_F_epsiWe  R_F_epsiHd R_F_epsiHk R_F_epsiEC R_F_epsiECAB] ;
% irf1(:,:,5)=[R_D_epsiA R_D_epsiJ R_D_epsiSm  R_D_epsiWe  R_D_epsiHd R_D_epsiHk R_D_epsiEC R_D_epsiECAB] ;
% irf1(:,:,6)=[Y_net_2_epsiA Y_net_2_epsiJ Y_net_2_epsiSm  Y_net_2_epsiWe  Y_net_2_epsiHd Y_net_2_epsiHk Y_net_2_epsiEC Y_net_2_epsiECAB] ;
% irf1(:,:,7)=[C_epsiA C_epsiJ C_epsiSm   C_epsiWe  C_epsiHd C_epsiHk C_epsiEC C_epsiECAB] ;
% irf1(:,:,8)=[I_epsiA I_epsiJ I_epsiSm   I_epsiWe  I_epsiHd I_epsiHk I_epsiEC I_epsiECAB] ;
% irf1(:,:,9)=[q_H_epsiA q_H_epsiJ q_H_epsiSm   q_H_epsiWe  q_H_epsiHd q_H_epsiHk q_H_epsiEC q_H_epsiECAB] ;
% 

irf1(:,:,1)=[b_m_obs_epsiA b_m_obs_epsiJ b_m_obs_epsiEC b_m_obs_epsiECAB b_m_obs_epsiLTVH] ;
irf1(:,:,2)=[b_e_obs_epsiA b_e_obs_epsiJ   b_e_obs_epsiEC b_e_obs_epsiECAB b_e_obs_epsiLTVH] ;
irf1(:,:,3)=[R_m_obs_epsiA R_m_obs_epsiJ     R_m_obs_epsiEC R_m_obs_epsiECAB R_m_obs_epsiLTVH] ;
irf1(:,:,4)=[R_F_obs_epsiA R_F_obs_epsiJ     R_F_obs_epsiEC R_F_obs_epsiECAB R_F_obs_epsiLTVH] ;
irf1(:,:,5)=[Y_net_2_obs_epsiA Y_net_2_obs_epsiJ   Y_net_2_obs_epsiEC Y_net_2_obs_epsiECAB Y_net_2_obs_epsiLTVH] ;
irf1(:,:,6)=[C_obs_epsiA C_obs_epsiJ   C_obs_epsiEC C_obs_epsiECAB C_obs_epsiLTVH] ;
irf1(:,:,7)=[q_H_obs_epsiA q_H_obs_epsiJ    q_H_obs_epsiEC q_H_obs_epsiECAB q_H_obs_epsiECAB] ;
%irf1(:,:,8)=[w_obs_epsiA w_obs_epsiJ w_obs_epsiSm   w_obs_epsiWe  w_obs_epsiHd w_obs_epsiHk w_obs_epsiEC w_obs_epsiECAB] ;



shock_names={'A','J','Se','We','Hd','Hk','EC','CAB'};

save irf1.mat irf1;

%  set_zeta_m=  0.6383 ; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
% set_zeta_F= 0.7820; %SS*** Interest rate stickiness,  appears in 51 FOC for bank business lending, 51 is not in use  %60
%  

zeta_m=0.15; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
zeta_F= 0.4121; %SS*** Interest rate stickiness,  appears in 51 FOC for bank business lending, 51 is not in use  %60

set_parameter_values_zeta(zeta_m,zeta_F);
dynare LTV1.mod%onlyclearglobals nolog nograph nointeractive;


% irf2(:,:,1)=[b_m_epsiA b_m_epsiJ b_m_epsiSm  b_m_epsiWe  b_m_epsiHd b_m_epsiHk b_m_epsiEC b_m_epsiECAB] ;
% irf2(:,:,2)=[b_e_epsiA b_e_epsiJ b_e_epsiSm   b_e_epsiWe  b_e_epsiHd b_e_epsiHk b_e_epsiEC b_e_epsiECAB] ;
% irf2(:,:,3)=[R_m_epsiA R_m_epsiJ R_m_epsiSm   R_m_epsiWe  R_m_epsiHd R_m_epsiHk R_m_epsiEC R_m_epsiECAB] ;
% irf2(:,:,4)=[R_F_epsiA R_F_epsiJ R_F_epsiSm   R_F_epsiWe  R_F_epsiHd R_F_epsiHk R_F_epsiEC R_F_epsiECAB] ;
% irf2(:,:,5)=[R_D_epsiA R_D_epsiJ R_D_epsiSm   R_D_epsiWe  R_D_epsiHd R_D_epsiHk R_D_epsiEC R_D_epsiECAB] ;
% irf2(:,:,6)=[Y_net_2_epsiA Y_net_2_epsiJ Y_net_2_epsiSm   Y_net_2_epsiWe  Y_net_2_epsiHd Y_net_2_epsiHk Y_net_2_epsiEC Y_net_2_epsiECAB] ;
% irf2(:,:,7)=[C_epsiA C_epsiJ C_epsiSm   C_epsiWe  C_epsiHd C_epsiHk C_epsiEC C_epsiECAB] ;
% irf2(:,:,8)=[I_epsiA I_epsiJ I_epsiSm   I_epsiWe  I_epsiHd I_epsiHk I_epsiEC I_epsiECAB] ;
% irf2(:,:,9)=[q_H_epsiA q_H_epsiJ q_H_epsiSm   q_H_epsiWe  q_H_epsiHd q_H_epsiHk q_H_epsiEC q_H_epsiECAB] ;

irf2(:,:,1)=[b_m_obs_epsiA b_m_obs_epsiJ b_m_obs_epsiEC b_m_obs_epsiECAB b_m_obs_epsiLTVH] ;
irf2(:,:,2)=[b_e_obs_epsiA b_e_obs_epsiJ   b_e_obs_epsiEC b_e_obs_epsiECAB b_e_obs_epsiLTVH] ;
irf2(:,:,3)=[R_m_obs_epsiA R_m_obs_epsiJ     R_m_obs_epsiEC R_m_obs_epsiECAB R_m_obs_epsiLTVH] ;
irf2(:,:,4)=[R_F_obs_epsiA R_F_obs_epsiJ     R_F_obs_epsiEC R_F_obs_epsiECAB R_F_obs_epsiLTVH] ;
irf2(:,:,5)=[Y_net_2_obs_epsiA Y_net_2_obs_epsiJ   Y_net_2_obs_epsiEC Y_net_2_obs_epsiECAB Y_net_2_obs_epsiLTVH] ;
irf2(:,:,6)=[C_obs_epsiA C_obs_epsiJ   C_obs_epsiEC C_obs_epsiECAB C_obs_epsiLTVH] ;
irf2(:,:,7)=[q_H_obs_epsiA q_H_obs_epsiJ    q_H_obs_epsiEC q_H_obs_epsiECAB q_H_obs_epsiECAB] ;



save irf2.mat irf2;

%  set_zeta_m=  0.6383 ; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
% set_zeta_F= 0.7820; %SS*** Interest rate stickiness,  appears in 51 FOC for bank business lending, 51 is not in use  %60
%  
zeta_m=0.6720; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
zeta_F= 0.15; %SS*** Interest rate stickiness,  appears in 51 FOC for bank business lending, 51 is not in use  %60

set_parameter_values_zeta(zeta_m,zeta_F);
dynare LTV1.mod%onlyclearglobals nolog nograph nointeractive;

% irf3(:,:,1)=[b_m_epsiA b_m_epsiJ b_m_epsiSm  b_m_epsiWe  b_m_epsiHd b_m_epsiHk b_m_epsiEC b_m_epsiECAB] ;
% irf3(:,:,2)=[b_e_epsiA b_e_epsiJ b_e_epsiSm   b_e_epsiWe  b_e_epsiHd b_e_epsiHk b_e_epsiEC b_e_epsiECAB] ;
% irf3(:,:,3)=[R_m_epsiA R_m_epsiJ R_m_epsiSm   R_m_epsiWe  R_m_epsiHd R_m_epsiHk R_m_epsiEC R_m_epsiECAB] ;
% irf3(:,:,4)=[R_F_epsiA R_F_epsiJ R_F_epsiSm   R_F_epsiWe  R_F_epsiHd R_F_epsiHk R_F_epsiEC R_F_epsiECAB] ;
% irf3(:,:,5)=[R_D_epsiA R_D_epsiJ R_D_epsiSm   R_D_epsiWe  R_D_epsiHd R_D_epsiHk R_D_epsiEC R_D_epsiECAB] ;
% irf3(:,:,6)=[Y_net_2_epsiA Y_net_2_epsiJ Y_net_2_epsiSm   Y_net_2_epsiWe  Y_net_2_epsiHd Y_net_2_epsiHk Y_net_2_epsiEC Y_net_2_epsiECAB] ;
% irf3(:,:,7)=[C_epsiA C_epsiJ C_epsiSm   C_epsiWe  C_epsiHd C_epsiHk C_epsiEC C_epsiECAB] ;
% irf3(:,:,8)=[I_epsiA I_epsiJ I_epsiSm   I_epsiWe  I_epsiHd I_epsiHk I_epsiEC I_epsiECAB] ;
% irf3(:,:,9)=[q_H_epsiA q_H_epsiJ q_H_epsiSm   q_H_epsiWe  q_H_epsiHd q_H_epsiHk q_H_epsiEC q_H_epsiECAB] ;

irf3(:,:,1)=[b_m_obs_epsiA b_m_obs_epsiJ b_m_obs_epsiEC b_m_obs_epsiECAB b_m_obs_epsiLTVH] ;
irf3(:,:,2)=[b_e_obs_epsiA b_e_obs_epsiJ   b_e_obs_epsiEC b_e_obs_epsiECAB b_e_obs_epsiLTVH] ;
irf3(:,:,3)=[R_m_obs_epsiA R_m_obs_epsiJ     R_m_obs_epsiEC R_m_obs_epsiECAB R_m_obs_epsiLTVH] ;
irf3(:,:,4)=[R_F_obs_epsiA R_F_obs_epsiJ     R_F_obs_epsiEC R_F_obs_epsiECAB R_F_obs_epsiLTVH] ;
irf3(:,:,5)=[Y_net_2_obs_epsiA Y_net_2_obs_epsiJ   Y_net_2_obs_epsiEC Y_net_2_obs_epsiECAB Y_net_2_obs_epsiLTVH] ;
irf3(:,:,6)=[C_obs_epsiA C_obs_epsiJ   C_obs_epsiEC C_obs_epsiECAB C_obs_epsiLTVH] ;
irf3(:,:,7)=[q_H_obs_epsiA q_H_obs_epsiJ    q_H_obs_epsiEC q_H_obs_epsiECAB q_H_obs_epsiECAB] ;
%irf3(:,:,8)=[w_obs_epsiA w_obs_epsiJ w_obs_epsiSm   w_obs_epsiWe  w_obs_epsiHd w_obs_epsiHk w_obs_epsiEC w_obs_epsiECAB] ;



save irf3.mat irf3;
clear;
%====================================================================\
zeta_m=0.15; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
zeta_F= 0.15; %SS*** Interest rate stickiness,  appears in 51 FOC for bank business lending, 51 is not in use  %60

set_parameter_values_zeta(zeta_m,zeta_F);
dynare LTV1.mod%onlyclearglobals nolog nograph nointeractive;

% irf4(:,:,1)=[b_m_epsiA b_m_epsiJ b_m_epsiSm  b_m_epsiWe  b_m_epsiHd b_m_epsiHk b_m_epsiEC b_m_epsiECAB] ;
% irf4(:,:,2)=[b_e_epsiA b_e_epsiJ b_e_epsiSm   b_e_epsiWe  b_e_epsiHd b_e_epsiHk b_e_epsiEC b_e_epsiECAB] ;
% irf4(:,:,3)=[R_m_epsiA R_m_epsiJ R_m_epsiSm   R_m_epsiWe  R_m_epsiHd R_m_epsiHk R_m_epsiEC R_m_epsiECAB] ;
% irf4(:,:,4)=[R_F_epsiA R_F_epsiJ R_F_epsiSm   R_F_epsiWe  R_F_epsiHd R_F_epsiHk R_F_epsiEC R_F_epsiECAB] ;
% irf4(:,:,5)=[R_D_epsiA R_D_epsiJ R_D_epsiSm   R_D_epsiWe  R_D_epsiHd R_D_epsiHk R_D_epsiEC R_D_epsiECAB] ;
% irf4(:,:,6)=[Y_net_2_epsiA Y_net_2_epsiJ Y_net_2_epsiSm   Y_net_2_epsiWe  Y_net_2_epsiHd Y_net_2_epsiHk Y_net_2_epsiEC Y_net_2_epsiECAB] ;
% irf4(:,:,7)=[C_epsiA C_epsiJ C_epsiSm   C_epsiWe  C_epsiHd C_epsiHk C_epsiEC C_epsiECAB] ;
% irf4(:,:,8)=[I_epsiA I_epsiJ I_epsiSm   I_epsiWe  I_epsiHd I_epsiHk I_epsiEC I_epsiECAB] ;
% irf4(:,:,9)=[q_H_epsiA q_H_epsiJ q_H_epsiSm   q_H_epsiWe  q_H_epsiHd q_H_epsiHk q_H_epsiEC q_H_epsiECAB] ;

irf4(:,:,1)=[b_m_obs_epsiA b_m_obs_epsiJ b_m_obs_epsiEC b_m_obs_epsiECAB b_m_obs_epsiLTVH] ;
irf4(:,:,2)=[b_e_obs_epsiA b_e_obs_epsiJ   b_e_obs_epsiEC b_e_obs_epsiECAB b_e_obs_epsiLTVH] ;
irf4(:,:,3)=[R_m_obs_epsiA R_m_obs_epsiJ     R_m_obs_epsiEC R_m_obs_epsiECAB R_m_obs_epsiLTVH] ;
irf4(:,:,4)=[R_F_obs_epsiA R_F_obs_epsiJ     R_F_obs_epsiEC R_F_obs_epsiECAB R_F_obs_epsiLTVH] ;
irf4(:,:,5)=[Y_net_2_obs_epsiA Y_net_2_obs_epsiJ   Y_net_2_obs_epsiEC Y_net_2_obs_epsiECAB Y_net_2_obs_epsiLTVH] ;
irf4(:,:,6)=[C_obs_epsiA C_obs_epsiJ   C_obs_epsiEC C_obs_epsiECAB C_obs_epsiLTVH] ;
irf4(:,:,7)=[q_H_obs_epsiA q_H_obs_epsiJ    q_H_obs_epsiEC q_H_obs_epsiECAB q_H_obs_epsiECAB] ;




save irf4.mat irf4;
clear;
%====================================================================


load irf1.mat;
load irf2.mat;
load irf3.mat;
load irf4.mat;

% irf1=irf1/100;
% irf2=irf2/100;
% irf3=irf3/100;
% irf4=irf4/100;


var_names_plot={'Household Borrowing','Business Borrowing','Household Interest Rate',...
    'Business Interest Rate','Aggegate Output','Aggregate Consumption','House Prices'};
shock_names={'A','J','Sm','We','Hd','Hk','EC','CAB'};
% shock_names_plot={'TFP (A)','Housing Preference (J)', 'Household Uncertainty Shock (Sm)','Business Net Worth (We)',...
%  'Housing Depreciation Shock (Hd)',...
%     'Capital Depreciation Shock (Hk)' 'Consumption Preference Shock (C)', 'Bank Capital Shock (CAB)'};

shock_names_plot={'TFP (A)','Housing Preference (J)',  'Consumption Preference Shock (C)', 'Bank Capital Shock (CAB)' 'Household LTV (LTVH)'};


%irf order: irf_length, shock, variable
for vv=1:length(shock_names_plot);
    f=figure;
    f.Name=char(shock_names_plot(vv));
    f.Units='normalized';
    f.Position=[0 -0.2 1 1];
    for ss=1:length(var_names_plot);
        subplot(2,4,ss);
       
        plot(irf1(:,vv,ss),'color','black');
        hold on;
          plot(irf2(:,vv,ss),'.-','color','red');
          hold on;
          plot(irf3(:,vv,ss),'--','color','green');
          hold on;
          plot(irf4(:,vv,ss),'color','blue');
          
           title(var_names_plot(ss));
    end
%      set_zeta_m= [0.691226178352873]; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
% set_zeta_F=  [0.447236412060053]; 
    legend('\zeta_m=0.69, \zeta_e=0.45',...
       '\zeta_m=0.1, \zeta_e=0.45',...
        '\zeta_m=0.69, \zeta_e=0.1',...
        '\zeta_m=0.1, \zeta_e=0.1');
        fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,['stickiness' char(shock_names(vv))],'-dpdf');
end
