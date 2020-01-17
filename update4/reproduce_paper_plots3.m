%TURN OFF zeta_m and zeta_F in estimation block of .mod file before running

clear;clc;%close all;
 Cyphi_H=0;Cyphi_F=0;phi_Hs=0.11;phi_Fs=0.11;phis=0.11;epsilonH1s=0.86;epsilonF1s=0.86;LTVHrule=0;LTVFrule=0;
%Cyphi_H=0;Cyphi_F=0;phi_Hs=0.11;phi_Fs=0.11;phis=0.11;epsilonH1s=0.86;epsilonF1s=0.86;LTVHrule=0;LTVFrule=0;
 zeta_m=0.691226178352873; 
zeta_F=0.447236412060053;
set_parameter_values_counterfactual( Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,phis,epsilonH1s,epsilonF1s,LTVHrule,LTVFrule,zeta_m,zeta_F)

dynare LTV1.mod noclearall;
irf_length=40;
numVar=7;
numShock=9;
%irf=nan(irf_length,numShock,numVar);


% irf1(:,:,1)=[b_m_epsiA b_m_epsiJ b_m_epsiSB b_m_epsiHd b_m_epsiECAB];
% irf1(:,:,2)=[b_e_epsiA b_e_epsiJ b_e_epsiSB b_e_epsiHd b_e_epsiECAB];
% irf1(:,:,3)=[R_m_epsiA R_m_epsiJ R_m_epsiSB R_m_epsiHd R_m_epsiECAB];
% irf1(:,:,4)=[R_F_epsiA R_F_epsiJ R_F_epsiSB R_F_epsiHd R_F_epsiECAB];
% irf1(:,:,5)=[q_H_epsiA q_H_epsiJ q_H_epsiSB q_H_epsiHd q_H_epsiECAB];
% irf1(:,:,6)=[Y_net_2_epsiA Y_net_2_epsiJ Y_net_2_epsiSB Y_net_2_epsiHd Y_net_2_epsiECAB];

irf1(:,:,1)=[b_m_obs_epsiA b_m_obs_epsiJ b_m_obs_epsiSB b_m_obs_epsiHd b_m_obs_epsiECAB b_m_obs_epsiHk];
irf1(:,:,2)=[b_e_obs_epsiA b_e_obs_epsiJ b_e_obs_epsiSB b_e_obs_epsiHd b_e_obs_epsiECAB b_e_obs_epsiHk];
irf1(:,:,3)=[R_m_obs_epsiA R_m_obs_epsiJ R_m_obs_epsiSB R_m_obs_epsiHd R_m_obs_epsiECAB R_m_obs_epsiHk];
irf1(:,:,4)=[R_F_obs_epsiA R_F_obs_epsiJ R_F_obs_epsiSB R_F_obs_epsiHd R_F_obs_epsiECAB R_F_obs_epsiHk];
irf1(:,:,5)=[q_H_obs_epsiA q_H_obs_epsiJ q_H_obs_epsiSB q_H_obs_epsiHd q_H_obs_epsiECAB q_H_obs_epsiHk];
irf1(:,:,6)=[Y_net_2_obs_epsiA Y_net_2_obs_epsiJ Y_net_2_obs_epsiSB Y_net_2_obs_epsiHd Y_net_2_obs_epsiECAB Y_net_2_obs_epsiHk];


% shock_names={'A','J','SB','Hd'};

save irf1.mat irf1;
 Cyphi_H=0;Cyphi_F=0;phi_Hs=0.15;phi_Fs=0.11;phis=0.11;epsilonH1s=0.86;epsilonF1s=0.86;LTVHrule=0;LTVFrule=0;
% Cyphi_H=0;Cyphi_F=0;phi_Hs=0.11;0;phi_Fs=0.11;0;phis=0.11;epsilonH1s=0.9;epsilonF1s=0.86;LTVHrule=0;LTVFrule=0;
%Cyphi_H=0.4;Cyphi_F=0;phi_Hs=0.11;0;phi_Fs=0.11;0;phis=0.11;epsilonH1s=0.86;epsilonF1s=0.86;LTVHrule=0;LTVFrule=0;
  

zeta_m=0.691226178352873; 
zeta_F=0.447236412060053;

set_parameter_values_counterfactual( Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,phis,epsilonH1s,epsilonF1s,LTVHrule,LTVFrule,zeta_m,zeta_F)
dynare LTV1.mod noclearall;%onlyclearglobals nolog nograph nointeractive;

% irf2(:,:,1)=[b_m_epsiA b_m_epsiJ b_m_epsiSB b_m_epsiHd b_m_epsiECAB];
% irf2(:,:,2)=[b_e_epsiA b_e_epsiJ b_e_epsiSB b_e_epsiHd b_e_epsiECAB];
% irf2(:,:,3)=[R_m_epsiA R_m_epsiJ R_m_epsiSB R_m_epsiHd R_m_epsiECAB];
% irf2(:,:,4)=[R_F_epsiA R_F_epsiJ R_F_epsiSB R_F_epsiHd R_F_epsiECAB];
% irf2(:,:,5)=[q_H_epsiA q_H_epsiJ q_H_epsiSB q_H_epsiHd q_H_epsiECAB];
% irf2(:,:,6)=[Y_net_2_epsiA Y_net_2_epsiJ Y_net_2_epsiSB Y_net_2_epsiHd Y_net_2_epsiECAB];


irf2(:,:,1)=[b_m_obs_epsiA b_m_obs_epsiJ b_m_obs_epsiSB b_m_obs_epsiHd b_m_obs_epsiECAB b_m_obs_epsiHk];
irf2(:,:,2)=[b_e_obs_epsiA b_e_obs_epsiJ b_e_obs_epsiSB b_e_obs_epsiHd b_e_obs_epsiECAB b_e_obs_epsiHk];
irf2(:,:,3)=[R_m_obs_epsiA R_m_obs_epsiJ R_m_obs_epsiSB R_m_obs_epsiHd R_m_obs_epsiECAB R_m_obs_epsiHk];
irf2(:,:,4)=[R_F_obs_epsiA R_F_obs_epsiJ R_F_obs_epsiSB R_F_obs_epsiHd R_F_obs_epsiECAB R_F_obs_epsiHk];
irf2(:,:,5)=[q_H_obs_epsiA q_H_obs_epsiJ q_H_obs_epsiSB q_H_obs_epsiHd q_H_obs_epsiECAB q_H_obs_epsiHk];
irf2(:,:,6)=[Y_net_2_obs_epsiA Y_net_2_obs_epsiJ Y_net_2_obs_epsiSB Y_net_2_obs_epsiHd Y_net_2_obs_epsiECAB Y_net_2_obs_epsiHk];

save irf2.mat irf2;


Cyphi_H=0;Cyphi_F=0;phi_Hs=0.11;phi_Fs=0.11;phis=0.11;epsilonH1s=0.86;epsilonF1s=0.86;LTVHrule=0;LTVFrule=0;
 zeta_m=0.1; 
zeta_F=0.1;

set_parameter_values_counterfactual( Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,phis,epsilonH1s,epsilonF1s,LTVHrule,LTVFrule,zeta_m,zeta_F)
dynare LTV1.mod noclearall;


% irf3(:,:,1)=[b_m_epsiA b_m_epsiJ b_m_epsiSB b_m_epsiHd b_m_epsiECAB];
% irf3(:,:,2)=[b_e_epsiA b_e_epsiJ b_e_epsiSB b_e_epsiHd b_e_epsiECAB];
% irf3(:,:,3)=[R_m_epsiA R_m_epsiJ R_m_epsiSB R_m_epsiHd R_m_epsiECAB];
% irf3(:,:,4)=[R_F_epsiA R_F_epsiJ R_F_epsiSB R_F_epsiHd R_F_epsiECAB];
% irf3(:,:,5)=[q_H_epsiA q_H_epsiJ q_H_epsiSB q_H_epsiHd q_H_epsiECAB];
% irf3(:,:,6)=[Y_net_2_epsiA Y_net_2_epsiJ Y_net_2_epsiSB Y_net_2_epsiHd Y_net_2_epsiECAB];

irf3(:,:,1)=[b_m_obs_epsiA b_m_obs_epsiJ b_m_obs_epsiSB b_m_obs_epsiHd b_m_obs_epsiECAB b_m_obs_epsiHk];
irf3(:,:,2)=[b_e_obs_epsiA b_e_obs_epsiJ b_e_obs_epsiSB b_e_obs_epsiHd b_e_obs_epsiECAB b_e_obs_epsiHk];
irf3(:,:,3)=[R_m_obs_epsiA R_m_obs_epsiJ R_m_obs_epsiSB R_m_obs_epsiHd R_m_obs_epsiECAB R_m_obs_epsiHk];
irf3(:,:,4)=[R_F_obs_epsiA R_F_obs_epsiJ R_F_obs_epsiSB R_F_obs_epsiHd R_F_obs_epsiECAB R_F_obs_epsiHk];
irf3(:,:,5)=[q_H_obs_epsiA q_H_obs_epsiJ q_H_obs_epsiSB q_H_obs_epsiHd q_H_obs_epsiECAB q_H_obs_epsiHk];
irf3(:,:,6)=[Y_net_2_obs_epsiA Y_net_2_obs_epsiJ Y_net_2_obs_epsiSB Y_net_2_obs_epsiHd Y_net_2_obs_epsiECAB Y_net_2_obs_epsiHk];

save irf3.mat irf3;

%================================
%Cyphi_H=0;Cyphi_F=0;phi_Hs=0.11;phi_Fs=0.11;phis=0;epsilonH1s=0.9;epsilonF1s=0.86;LTVHrule=0;LTVFrule=0;
% Cyphi_H=0.4;Cyphi_F=0.4;phi_Hs=0.11;phi_Fs=0.11;phis=0.11;epsilonH1s=0.86;epsilonF1s=0.86;LTVHrule=0;LTVFrule=0;
 Cyphi_H=0;Cyphi_F=0;phi_Hs=0.15;0;phi_Fs=0.11;0;phis=0;epsilonH1s=0.86;epsilonF1s=0.86;LTVHrule=0;LTVFrule=0;
 
zeta_m=0.1; 
zeta_F=0.1;

set_parameter_values_counterfactual( Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,phis,epsilonH1s,epsilonF1s,LTVHrule,LTVFrule,zeta_m,zeta_F)
dynare LTV1.mod noclearall;


% irf4(:,:,1)=[b_m_epsiA b_m_epsiJ b_m_epsiSB b_m_epsiHd b_m_epsiECAB];
% irf4(:,:,2)=[b_e_epsiA b_e_epsiJ b_e_epsiSB b_e_epsiHd b_e_epsiECAB];
% irf4(:,:,3)=[R_m_epsiA R_m_epsiJ R_m_epsiSB R_m_epsiHd R_m_epsiECAB];
% irf4(:,:,4)=[R_F_epsiA R_F_epsiJ R_F_epsiSB R_F_epsiHd R_F_epsiECAB];
% irf4(:,:,5)=[q_H_epsiA q_H_epsiJ q_H_epsiSB q_H_epsiHd q_H_epsiECAB];
% irf4(:,:,6)=[Y_net_2_epsiA Y_net_2_epsiJ Y_net_2_epsiSB Y_net_2_epsiHd Y_net_2_epsiECAB];

irf4(:,:,1)=[b_m_obs_epsiA b_m_obs_epsiJ b_m_obs_epsiSB b_m_obs_epsiHd b_m_obs_epsiECAB b_m_obs_epsiHk];
irf4(:,:,2)=[b_e_obs_epsiA b_e_obs_epsiJ b_e_obs_epsiSB b_e_obs_epsiHd b_e_obs_epsiECAB b_e_obs_epsiHk];
irf4(:,:,3)=[R_m_obs_epsiA R_m_obs_epsiJ R_m_obs_epsiSB R_m_obs_epsiHd R_m_obs_epsiECAB R_m_obs_epsiHk];
irf4(:,:,4)=[R_F_obs_epsiA R_F_obs_epsiJ R_F_obs_epsiSB R_F_obs_epsiHd R_F_obs_epsiECAB R_F_obs_epsiHk];
irf4(:,:,5)=[q_H_obs_epsiA q_H_obs_epsiJ q_H_obs_epsiSB q_H_obs_epsiHd q_H_obs_epsiECAB q_H_obs_epsiHk];
irf4(:,:,6)=[Y_net_2_obs_epsiA Y_net_2_obs_epsiJ Y_net_2_obs_epsiSB Y_net_2_obs_epsiHd Y_net_2_obs_epsiECAB Y_net_2_obs_epsiHk];


save irf4.mat irf4;



%====================================================================
clear;

load irf1.mat;
load irf2.mat;
load irf3.mat;
load irf4.mat;

% irf1=irf1/100;
% irf2=irf2/100;
% irf3=irf3/100;
% irf4=irf4/100;

var_names={'b_m' 'b_e' 'R_m' 'R_F' 'q_H' 'Y_net_2'};
var_names_plot={'Mortgage Lending','Corporate Lending','Household Interest Rate',...
    'Business Interest Rate','House Price','Aggregate Output'};
shock_names={'A','J','SB','HD','ECAB', 'Hk'};
shock_names_plot={'TFP (A)','Housing Preference (J)', 'Bank Uncertainty Shock (SB)',...
 'Housing Depreciation Shock (Hd)', 'Bank Capital Shock (ECAB)', 'Capital Depreciation Shock (Hk)'};


%irf order: irf_length, shock, variable
for vv=1:length(shock_names);
% for vv=1:2;
    f=figure;
    f.Name=char(shock_names_plot(vv));
    f.Units='normalized';
    f.Position=[0.2 0.2 1 0.5];
%     for ss=1:length(var_names);
for ss=1:2
          subplot(1,2,ss);      

plot(irf1(:,vv,ss),'color','black');
        hold on;
          plot(irf2(:,vv,ss),'.-','color','red');
          hold on;
          plot(irf3(:,vv,ss),'--','color','green');
          hold on;
          plot(irf4(:,vv,ss),'color','blue');
          
           title(var_names_plot(ss));
           
                   if ss==1
%             legend('Estimated Stickiness, No CCyB',...
%         'Estimated Stickiness, with CCyB',...
%         'Low Stickiness, No CcyB',...
%         'Low Stickiness, with CCyB');

%             legend('Estimated Stickiness, LTV=0.86',...
%         'Estimated Stickiness, LTV=0.9',...
%         'Low Stickiness, LTV=0.86',...
%         'Low Stickiness, LTV=0.9');

            legend('Estimated Stickiness, Low CAR',...
        'Estimated Stickiness, High CAR',...
        'Low Stickiness, Low CAR',...
        'Low Stickiness, High CAR');

        end
           
    end
%      set_zeta_m= [0.691226178352873]; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
% set_zeta_F=  [0.447236412060053]; 

        fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,['stickiness_CAR' char(shock_names(vv))],'-dpdf');
end

%1&2 with stickiness 
%3&4 without stickiness 
%irf order: irf_length, shock, variable
%bm=1,be=2,hd=4,ecab=5

irf1_cumsum=cumsum(irf1);
irf2_cumsum=cumsum(irf2);
irf3_cumsum=cumsum(irf3);
irf4_cumsum=cumsum(irf4);

diff_sticky=cumsum(abs(irf1-irf2));
diff_nosticky=cumsum(abs(irf3-irf4));

% %irf order: irf_length, shock, variable
% for vv=1:length(shock_names);
% % for vv=1:2;
%     f=figure;
%     f.Name=char(shock_names_plot(vv));
%     f.Units='normalized';
%     f.Position=[0.2 0.2 1 0.5];
% %     for ss=1:length(var_names);
% for ss=1:2
%           subplot(1,2,ss);      
% 
% plot(irf1_cumsum(:,vv,ss),'color','black');
%         hold on;
%           plot(irf2_cumsum(:,vv,ss),'.-','color','red');
%           hold on;
%           plot(irf3_cumsum(:,vv,ss),'--','color','green');
%           hold on;
%           plot(irf4_cumsum(:,vv,ss),'color','blue');
%           
%            title(var_names_plot(ss));
%            
%                    if ss==1
% %             legend('Estimated Stickiness, No CCyB',...
% %         'Estimated Stickiness, with CCyB',...
% %         'Low Stickiness, No CcyB',...
% %         'Low Stickiness, with CCyB');
% 
% %             legend('Estimated Stickiness, LTV=0.86',...
% %         'Estimated Stickiness, LTV=0.9',...
% %         'Low Stickiness, LTV=0.86',...
% %         'Low Stickiness, LTV=0.9');
% 
%             legend('Estimated Stickiness, Low CAR',...
%         'Estimated Stickiness, High CAR',...
%         'Low Stickiness, Low CAR',...
%         'Low Stickiness, High CAR');
% 
%         end
%            
%     end
% %      set_zeta_m= [0.691226178352873]; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
% % set_zeta_F=  [0.447236412060053]; 
% 
%         fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,['stickiness_cumul_CAR' char(shock_names(vv))],'-dpdf');
% end

% diff_sticky=abs(irf1_cumsum-irf2_cumsum);
% diff_nosticky=abs(irf3_cumsum-irf4_cumsum);
b_m_J=[-diff_sticky(end,2,1)+diff_nosticky(end,2,1)]
b_e_J=[-diff_sticky(end,2,2)+diff_nosticky(end,2,2)]
b_m_HD=[-diff_sticky(end,4,1)+diff_nosticky(end,4,1)]
b_m_ECAB=[-diff_sticky(end,5,1)+diff_nosticky(end,5,1)]
b_m_Hk=[-diff_sticky(end,6,1)+diff_nosticky(end,6,1)]
b_e_HD=[-diff_sticky(end,4,2)+diff_nosticky(end,4,2)]
b_e_ECAB=[-diff_sticky(end,5,2)+diff_nosticky(end,5,2)]
b_e_Hk=[-diff_sticky(end,6,2)+diff_nosticky(end,6,2)]