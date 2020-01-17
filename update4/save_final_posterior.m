clear;clc;close all;
%load LTV1_mode.mat;

%posteriorMode=xparam1;
%1--table name, 2--latex name
parameter_names_long={'St Dev Hd','$\epsilon_{Hd}$';...
    'St Dev Hk','$\epsilon_{Hk}$';...
    'St Dev A','$\epsilon_{A}$';...
    'St Dev J','$\epsilon_{J}$';...
    'St Dev Se','$\epsilon_{Se}$';...
    'St Dev Sm','$\epsilon_{Sm}$';...
    'St Dev SB','$\epsilon_{SB}$';...
    'St Dev We','$\epsilon_{We}$';...
    'St Dev Wb','$\epsilon_{Wb}$';...
    'St Dev Markup_m','$\epsilon_{markup_m}$';...
'St Dev Markup_e','$\epsilon_{markup_e}$';...
'St Dev EC','$\epsilon_{EC}$';...
'St Dev ECAB','$\epsilon_{ECAB}$';...
'\sigma_e','$\sigma_e$';...
'sigma_m','$\sigma_m$';...
'sigma_B','$\sigma_B$';...
'\psi_i','$\psi_i$';...
'\psi_h','$\psi_h$';...
'\zeta_m','$\zeta_m$';...
'\zeta_e','$\zeta_e$';...
'\rho_{Hd}','$\rho_{Hd}$';...
'\rho_{Hk}','$\rho_{Hk}$';...
'\rho_A','$\rho_{A}$';...
'\rho_J','$\rho_{J}$';...
'\rho_{Se}','$\rho_{Se}$';...
'\rho_{Sm}','$\rho_{Sm}$';...
'\rho_{SB}','$\rho_{SB}$';...
'\rho_{We}','$\rho_{We}$';...
'\rho_{Wb}','$\rho_{Wb}$';...
'\rho_{Markup_m}','$\rho_{markup_m}$';...
'\rho_{Markup_e}','$\rho_{markup_e}$';...
'\rho_{EC}','$\rho_{EC}$';...
'\rho_{ECAB}','$\rho_{ECAB}$'};



%modeHess=hh;
N=100000;
mh_conf_sig=0.95;
round_dc=4;
%numVar=length(posteriorMode);
numVar=31;
%load each file and save the second half 
%load('c:\users\tolga\desktop\3d_model_tolga\version_oct_15\LTV1\metropolis\LTV1_mh1_blck1.mat');
load('c:\users\tolga\desktop\3d_model_updated\main_codes\LTV1\metropolis\LTV1_mh1_blck1.mat');

logpost_final(1:N)=logpo2(N+1:end);
posteriorDist(1:N,:)=x2(N+1:end,:);

load('c:\users\tolga\desktop\3d_model_updated\main_codes\LTV1\metropolis\LTV1_mh1_blck2.mat');
logpost_final(N+1:2*N)=logpo2(N+1:end);
posteriorDist(N+1:2*N,:)=x2(N+1:end,:);
% 
load('c:\users\tolga\desktop\3d_model_updated\main_codes\LTV1\metropolis\LTV1_mh1_blck3.mat');
logpost_final(2*N+1:3*N)=logpo2(N+1:end);
posteriorDist(2*N+1:3*N,:)=x2(N+1:end,:);
% 
load('c:\users\tolga\desktop\3d_model_updated\main_codes\LTV1\metropolis\LTV1_mh1_blck4.mat');
logpost_final(3*N+1:4*N)=logpo2(N+1:end);
posteriorDist(3*N+1:4*N,:)=x2(N+1:end,:);

save MH_draws.mat logpost_final posteriorDist;

posteriorMean=mean(posteriorDist)';
posteriorDist=sort(posteriorDist);
hpdDraws=((1-mh_conf_sig)*length(posteriorDist));


for paraInd=1:numVar
    jj=size(posteriorDist,1)-hpdDraws-2;

for ii=1:hpdDraws
    kk(ii,paraInd)=posteriorDist(jj,paraInd)-posteriorDist(ii,paraInd);
    jj=jj+1;
end
[kmin,idx]=min(kk(:,paraInd));
hpd_interval(paraInd,:)=[posteriorDist(idx,paraInd) posteriorDist(idx,paraInd)+kmin];
post_deciles(paraInd,:)=posteriorDist([(0.05*length(posteriorDist(:,paraInd)))...
    round(0.2*length(posteriorDist(:,paraInd)))...
    round(0.3*length(posteriorDist(:,paraInd)))...
    round(0.4*length(posteriorDist(:,paraInd)))...
    round(0.5*length(posteriorDist(:,paraInd)))...
    round(0.6*length(posteriorDist(:,paraInd)))...
    round(0.7*length(posteriorDist(:,paraInd)))...
    round(0.8*length(posteriorDist(:,paraInd)))...
    round(0.95*length(posteriorDist(:,paraInd)))],paraInd);
end


figure('Name','Posterior Distributions I','units','normalized','outerposition',[0 0 1 1]);
index=0;
for ii=1:16
    index=index+1;
%     plot_var=posteriorDist(:,ii)/sum(posteriorDist(:,ii));
plot_var=posteriorDist(:,ii);
    
subplot(3,6,index);hist(plot_var);
%hold on;
%plot([posteriorMode(ii) posteriorMode(ii)],[ylim],'color','blue'); 
hold on;
plot([posteriorMean(ii) posteriorMean(ii)],[ylim],'color','red'); 
title(parameter_names_long(ii,1))
end

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
string=['posteriordistributions1'];
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,string,'-dpdf');

figure('Name','Posterior Distributions II','units','normalized','outerposition',[0 0 1 1]);
index=0;
for ii=17:numVar
    index=index+1;
%     plot_var=posteriorDist(:,ii)/sum(posteriorDist(:,ii));
plot_var=posteriorDist(:,ii);
    
subplot(3,6,index);hist(plot_var);
hold on;
%plot([posteriorMode(ii) posteriorMode(ii)],[ylim],'color','blue'); 
%hold on;
plot([posteriorMean(ii) posteriorMean(ii)],[ylim],'color','red'); 
title(parameter_names_long(ii))
end

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
string=['posteriordistributions2'];
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,string,'-dpdf');



%plot the prior as well for these two
mu=0.5;sigma2=0.2^2;a=(1-mu)*mu*mu/sigma2-mu;b=a*(1/mu-1);
% priorzeta=@(x) betapdf(x,a,b);
priorzeta= betarnd(a,b,[length(posteriorDist(:,1)),1]);

figure('Name','Posterior Distributions Calvo','units','normalized','outerposition',[0 0 1 0.5]);
subplot(1,3,1);

% fplot(priorzeta,[0 1],'--','lineWidth',3,'color','red');
% fplot(priorzeta,[0 1],'color','black');
histogram(priorzeta,'Normalization','probability');
title('Prior Distribution');
index=1;
for ii=19:20
    index=index+1;
%     plot_var=posteriorDist(:,ii)/sum(posteriorDist(:,ii));
plot_var=posteriorDist(:,ii);
% [ff xx]=ksdensity(posteriorDist(:,ii))
    
subplot(1,3,index);
histogram(plot_var,'Normalization','probability');
% histogram(plot_var);
xlim([0 1]);
%  plot(xx,ff,'color','black','lineWidth',3);
% hold on;
% plot([posteriorMode(ii) posteriorMode(ii)],[ylim],'color','blue','lineWidth',3); 
% hold on;
% plot([posteriorMean(ii) posteriorMean(ii)],[ylim],'color','red'); 
% hold on;
% fplot(priorzeta,[0 1],'--','lineWidth',3,'color','red');
%  histogram(priorzeta,'Normalization','probability');
%title(parameter_names_long(ii))
if ii==19
title('Interest Rate Stickiness: Mortgage Sector');
else 
title('Interest Rate Stickiness: Corporate Sector');
end
% legend('Posterior','Prior');

end


fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
string=['posteriordistributions_calvo2'];
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,string,'-dpdf');



hpd_LB=post_deciles(:,1);
hpd_UB=post_deciles(:,end);

posteriorMode=round(posteriorMode,round_dc);
posteriorMean=round(posteriorMean,round_dc);
hpd_LB=round(hpd_LB,round_dc);
hpd_UB=round(hpd_UB,round_dc);

table(parameter_names_long(:,2),posteriorMode,posteriorMean, hpd_LB, hpd_UB)

  output_file='Estimation_Results.csv';
 output_sheet='Estimation_Results';
  xlswrite(output_file,parameter_names_long(:,2),output_sheet,'a3');
   xlswrite(output_file,posteriorMode,output_sheet,'b3');
 xlswrite(output_file,posteriorMean,output_sheet,'c3');
 xlswrite(output_file,hpd_LB,output_sheet,'d3');
    xlswrite(output_file,hpd_UB,output_sheet,'e3');
    



 