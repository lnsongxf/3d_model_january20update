clear;clc;%close all
% phi_Hs_grid=[0.11,0.35];
% epsilonH1s_grid=[0.95,0.7];
% epsilonH1s_grid=[0.86 0.76];
% phi_Hs_grid=[0.11,0.15];
% epsilonH1s_grid=[0.86,0.76];
epsilonH1s_grid=0.86;
% phi_Hs_grid=[0.11, 0.19];
phi_Hs_grid=0.15;
Cyphi_grid=[0,0.3];
%  legend('CAR-11%,LTV-86%, No CCyB',...
%      'CAR-11%,LTV-86%, CCyB=0.5',...
% 'CAR-11%,LTV-76%, No CCyB',...
%      'CAR-11%,LTV-76%, CCyB=0.5',...
%      'CAR-15%,LTV-86%, No CCyB',...
%       'CAR-15%,LTV-86%, CCyB=0.5',...
% 'CAR-15%,LTV-76%, No CCyB',...
%   'CAR-15%,LTV-76%, CCyB=0.5');
for jj=1:length(phi_Hs_grid)
    for ee=1:length(epsilonH1s_grid)
        for cc=1:length(Cyphi_grid)
phi_Fs=phi_Hs_grid(jj);
%phi_Hs=0.11;
    phi_Hs=phi_Hs_grid(jj);
%         phis=phi_Hs_grid(jj);
phis=0.11;
%  phi_Fs=0.11;
% phi_Hs=0.11;
        
%epsilonH1s=0.86;
     epsilonH1s=epsilonH1s_grid(ee);%0.86;
%epsilonF1s=epsilonH1s_grid(ee);%0.86;
 epsilonF1s=0.86;
% Cyphi_H=phi_Hs_grid(jj);
% % Cyphi_F=phi_Hs_grid(jj);
Cyphi_H=Cyphi_grid(cc);
% Cyphi_F=Cyphi_grid(cc);
Cyphi_F=0;


LTVHrule=0;
LTVFrule=0;
    
  
set_parameter_values_policy( Cyphi_H,Cyphi_F,phi_Hs,phi_Fs,phis,epsilonH1s,epsilonF1s,LTVHrule,LTVFrule);
    dynare LTV1.mod onlyclearglobals nolog nograph nointeractive;
                        
%    IRF(:,:,1,jj,ee)=[Y_net_epsiA,I_epsiA,IH_epsiA,def_rate_e_epsiA,def_rate_m_epsiA,C_epsiA,b_m_epsiA,b_e_epsiA,R_m_epsiA,R_D_epsiA]/100;                       %IRF(irf_length,variable,shock,phi_Hs_grid,epsilonH1s_grid)
%    IRF(:,:,2,jj,ee)=[Y_net_epsiHd,I_epsiHd,IH_epsiHd,def_rate_e_epsiHd,def_rate_m_epsiHd,C_epsiHd,b_m_epsiHd,b_e_epsiHd,R_m_epsiHd,R_D_epsiHd]/100;             
%    IRF(:,:,3,jj,ee)=[Y_net_epsiWb,I_epsiWb,IH_epsiWb,def_rate_e_epsiWb,def_rate_m_epsiWb,C_epsiWb,b_m_epsiWb,b_e_epsiWb,R_m_epsiWb,R_D_epsiWb]/100;             
%    IRF(:,:,4,jj,ee)=[Y_net_epsiSe,I_epsiSe,IH_epsiSe,def_rate_e_epsiSe,def_rate_m_epsiSe,C_epsiSe,b_m_epsiSe,b_e_epsiSe,R_m_epsiSe,R_D_epsiSe]/100;                     
%    IRF(:,:,5,jj,ee)=[Y_net_epsiHk,I_epsiHk,IH_epsiHk,def_rate_e_epsiHk,def_rate_m_epsiHk,C_epsiHk,b_m_epsiHk,b_e_epsiHk,R_m_epsiHk,R_D_epsiHk]/100;          

    IRF(:,:,1,jj,ee,cc)=[Y_net_obs_epsiA,I_obs_epsiA,b_m_obs_epsiA,b_e_obs_epsiA]/100;                       %IRF(irf_length,variable,shock,phi_Hs_grid,epsilonH1s_grid)
    IRF(:,:,2,jj,ee,cc)=[Y_net_obs_epsiHd,I_obs_epsiHd,b_m_obs_epsiHd,b_e_obs_epsiHd]/100;                  
    IRF(:,:,3,jj,ee,cc)=[Y_net_obs_epsiSm,I_obs_epsiSm,b_m_obs_epsiSm,b_e_obs_epsiSm]/100;                                
    IRF(:,:,4,jj,ee,cc)=[Y_net_obs_epsiHk,I_obs_epsiHk,b_m_obs_epsiHk,b_e_obs_epsiHk]/100;                
    IRF(:,:,5,jj,ee,cc)=[Y_net_obs_epsiECAB,I_obs_epsiECAB,b_m_obs_epsiECAB,b_e_obs_epsiECAB]/100;            
    IRF(:,:,6,jj,ee,cc)=[Y_net_obs_epsiJ,I_obs_epsiJ,b_m_obs_epsiJ,b_e_obs_epsiJ]/100;      

    end
    end
end

%============================================================
%============================================================
figure('Name','TFP Shock','units','normalized','outerposition',[0 0 1 1]);
irf_length=size(IRF,1);
numVar=size(IRF,2);
numShocks=size(IRF,3);
%var_names={'Output','Investment','Housing Investment','Def. Rate e','Def. Rate m','Consumption','HH debt','Business debt','HH interest rate','bank rate'};
var_names={'Output','Investment','Business Debt','Household Debt'};

for ii=1:numVar;
    subplot(2,2,ii);
for jj=1:length(phi_Hs_grid)
            for ee=1:length(epsilonH1s_grid)
                for cc=1:length(Cyphi_grid)
                    
   
        
                    plot(squeeze(IRF(:,ii,1,jj,ee,cc)),'lineWidth',3);
    
            hold on;
         end
            end
                end
if ii==numVar
 legend('CAR-11%,LTV-86%, No CCyB',...
     'CAR-11%,LTV-86%, CCyB=0.5',...
'CAR-11%,LTV-76%, No CCyB',...
     'CAR-11%,LTV-76%, CCyB=0.5',...
     'CAR-15%,LTV-86%, No CCyB',...
      'CAR-15%,LTV-86%, CCyB=0.5',...
'CAR-15%,LTV-76%, No CCyB',...
  'CAR-15%,LTV-76%, CCyB=0.5');
end
    title(var_names(ii));
end
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'irf_epsiA','-dpdf');

figure('Name','Housing Depreciation Shock','units','normalized','outerposition',[0 0 1 1]);
for ii=1:numVar;
    subplot(2,2,ii);
for jj=1:length(phi_Hs_grid)
            for ee=1:length(epsilonH1s_grid)
                for cc=1:length(Cyphi_grid)
                    
   
        
                    plot(squeeze(IRF(:,ii,2,jj,ee,cc)),'lineWidth',3);
    
            hold on;
         end
            end
                end
if ii==numVar
 legend('CAR-11%,LTV-86%, No CCyB',...
     'CAR-11%,LTV-86%, CCyB=0.5',...
'CAR-11%,LTV-76%, No CCyB',...
     'CAR-11%,LTV-76%, CCyB=0.5',...
     'CAR-15%,LTV-86%, No CCyB',...
      'CAR-15%,LTV-86%, CCyB=0.5',...
'CAR-15%,LTV-76%, No CCyB',...
  'CAR-15%,LTV-76%, CCyB=0.5');
end
    title(var_names(ii));
end
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'irf_epsiHd','-dpdf');



figure('Name','Household Risk Shock','units','normalized','outerposition',[0 0 1 1]);
for ii=1:numVar;
    subplot(2,2,ii);
for jj=1:length(phi_Hs_grid)
            for ee=1:length(epsilonH1s_grid)
                for cc=1:length(Cyphi_grid)
                    
   
        
                    plot(squeeze(IRF(:,ii,3,jj,ee,cc)),'lineWidth',3);
    
            hold on;
         end
            end
                end
if ii==numVar
 legend('CAR-11%,LTV-86%, No CCyB',...
     'CAR-11%,LTV-86%, CCyB=0.5',...
'CAR-11%,LTV-76%, No CCyB',...
     'CAR-11%,LTV-76%, CCyB=0.5',...
     'CAR-15%,LTV-86%, No CCyB',...
      'CAR-15%,LTV-86%, CCyB=0.5',...
'CAR-15%,LTV-76%, No CCyB',...
  'CAR-15%,LTV-76%, CCyB=0.5');
end
    title(var_names(ii));
end
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'irf_epsiSm','-dpdf');


figure('Name','Capital Depreciation Shock','units','normalized','outerposition',[0 0 1 1]);
for ii=1:numVar;
    subplot(2,2,ii);
for jj=1:length(phi_Hs_grid)
            for ee=1:length(epsilonH1s_grid)
                for cc=1:length(Cyphi_grid)
                    
   
        
                    plot(squeeze(IRF(:,ii,4,jj,ee,cc)),'lineWidth',3);
    
            hold on;
         end
            end
                end
if ii==numVar
 legend('CAR-11%,LTV-86%, No CCyB',...
     'CAR-11%,LTV-86%, CCyB=0.5',...
'CAR-11%,LTV-76%, No CCyB',...
     'CAR-11%,LTV-76%, CCyB=0.5',...
     'CAR-15%,LTV-86%, No CCyB',...
      'CAR-15%,LTV-86%, CCyB=0.5',...
'CAR-15%,LTV-76%, No CCyB',...
  'CAR-15%,LTV-76%, CCyB=0.5');
end
    title(var_names(ii));
end
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'irf_epsiHk','-dpdf');


figure('Name','Bank Capital Shock','units','normalized','outerposition',[0 0 1 1]);
for ii=1:numVar;
    subplot(2,2,ii);
for jj=1:length(phi_Hs_grid)
            for ee=1:length(epsilonH1s_grid)
                for cc=1:length(Cyphi_grid)
                    
   
        
                    plot(squeeze(IRF(:,ii,5,jj,ee,cc)),'lineWidth',3);
    
            hold on;
         end
            end
                end
if ii==numVar
 legend('CAR-11%,LTV-86%, No CCyB',...
     'CAR-11%,LTV-86%, CCyB=0.5',...
'CAR-11%,LTV-76%, No CCyB',...
     'CAR-11%,LTV-76%, CCyB=0.5',...
     'CAR-15%,LTV-86%, No CCyB',...
      'CAR-15%,LTV-86%, CCyB=0.5',...
'CAR-15%,LTV-76%, No CCyB',...
  'CAR-15%,LTV-76%, CCyB=0.5');
end
    title(var_names(ii));
end
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'irf_epsiECAB','-dpdf');


figure('Name','Housing Preference Shock','units','normalized','outerposition',[0 0 1 1]);
for ii=1:numVar;
    subplot(2,2,ii);
for jj=1:length(phi_Hs_grid)
            for ee=1:length(epsilonH1s_grid)
                for cc=1:length(Cyphi_grid)
                    
   
        
                    plot(squeeze(IRF(:,ii,6,jj,ee,cc)),'lineWidth',3);
    
            hold on;
         end
            end
                end
if ii==numVar
 legend('CAR-11%,LTV-86%, No CCyB',...
     'CAR-11%,LTV-86%, CCyB=0.5',...
'CAR-11%,LTV-76%, No CCyB',...
     'CAR-11%,LTV-76%, CCyB=0.5',...
     'CAR-15%,LTV-86%, No CCyB',...
      'CAR-15%,LTV-86%, CCyB=0.5',...
'CAR-15%,LTV-76%, No CCyB',...
  'CAR-15%,LTV-76%, CCyB=0.5');
end
    title(var_names(ii));
end
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'irf_epsiJ','-dpdf');
