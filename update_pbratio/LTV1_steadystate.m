
function [ys,check,input] = ThreeDAdj_RD_steadystate(ys,exo)
format long
check=0;

%%%%%**********************************************************************   
global M_ 

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% 
% load LTV1_parameter_values.mat;
% sigma_epsiHd=set_sigma_epsiHd;
% sigma_epsiHk=set_sigma_epsiHk;
%  sigma_epsiA= set_sigma_epsiA; 
%  sigma_epsiJ= set_sigma_epsiJ; 
%  sigma_epsiK = set_sigma_epsiK;
% sigma_epsiH =set_sigma_epsiH ; 
% sigma_epsiSe=set_sigma_epsiSe; 
% sigma_epsiSm=set_sigma_epsiSm; 
%   sigma_epsiSF=  set_sigma_epsiSF;
%  sigma_epsiSH= set_sigma_epsiSH; 
%   sigma_epsiWb=set_sigma_epsiWb;
%   sigma_epsiRW=set_sigma_epsiRW;
%   sigma_epsiWe=  set_sigma_epsiWe;
%  pp  = set_pp ;
% betta_s=set_betta_s; 
% betta_m=set_betta_m; 
% hab =set_hab;
% Cyphi_H =set_Cyphi_H ; 
% phi_Hs=set_phi_Hs; 
% Cyphi_F=set_Cyphi_F;
% phi_Fs =set_phi_Fs  ;
% alphaa=set_alphaa; 
% rp=set_rp;
% rpe=set_rpe; 
% delta_K= set_delta_K;
% delta_H=set_delta_H; 
% mu_m=set_mu_m;
% mu_e=set_mu_e;
% mu_B=set_mu_B;
% mu_B=set_mu_B;
% sigma_e1=set_sigma_e1;
% sigma_m1=set_sigma_m1;
% sigma_B=set_sigma_B;
% sigma_B=set_sigma_B;
% varphi_s=set_varphi_s;
%  varphi_m= set_varphi_m; 
% v_s=set_v_s; 
% v_m=set_v_m; 
% chi_b=set_chi_b; 
% chi_e=set_chi_e; 
%   eta = set_eta;
%  a_e =set_a_e;
% a_s=set_a_s; 
%  a_b = set_a_b; 
% psi_i=set_psi_i; 
% psi_h=set_psi_h;
% rhoA =set_rhoA;
% rhoJ =set_rhoJ ;
% rhoK =set_rhoK;
%  rhoSe = set_rhoSe;
% rhoSm =set_rhoSm; 
% rhoSF =set_rhoSF; 
% rhoSH =set_rhoSH; 
%  rhoWb = set_rhoWb; 
% rhoWe =set_rhoWe; 
% rhoHd =set_rhoHd; 
% rhoHk =set_rhoHk ; 
% rhoRW =set_rhoRW ; 
% zeta1=set_zeta1; 
% zetae=set_zetae;
% epsilonH1s=set_epsilonH1s;
% epsilonF1s=set_epsilonF1s;
% tau=set_tau; 
% taue=set_taue;
% phiinf1=set_phiinf1; 
%  kappa = set_kappa  ;
% nu=set_nu; 
%  psib = set_psib; 
% phis=set_phis; 
% LTVHrule=set_LTVHrule;
% LTVFrule=set_LTVFrule;
%  gamma_y=set_gamma_y;
%  gamma_db=set_gamma_db;
%%%%%**********************************************************************
%ERW=1;
phi_F=phi_Fs;
phi_H=phi_Hs;

phi=phis;
%ELTV=1;
epsilonH=epsilonH1s;
epsilonF=epsilonF1s;
R_DDs = (1/betta_s);%where is this coming from?

%rhos = 1/(1-chi_b);
%rho_Hs = rho_Fs;
R_Hs = 1-delta_H;
As=1;
% Capital adjustment costs = 0
q_Ks=1;
g_Is= 0;
g_I_1s=0 ;
PIs=0;
PHs=0;
g_Hs= 0;
g_H_1s=0 ;

%%%%%**********************************************************************
% 2. The remaining variables require the use of a solver (system of equations).
% However, by assuming that the steady state value of the following variables is known
% R_Ks = xfs(1);       % Rate of return to capital
% Ls   = xfs(2);       % Aggregate hours worked
% L_ms = xfs(3);       % Labour of borrowers
% x_es = xfs(4);       % Corporate leverage
% x_ms = xfs(5);       % Household leverage        		 
% omega_bar_Bs=xfs(6); % Mortgage bank default cut off 
% R_tilde_Hs=xfs(7);   % Return on portfolio of mortgage loans
% omega_bar_Bs=xfs(8); % Corporate bank default cut off
% R_tilde_Fs=xfs(9);   % Return on portfolio of corporate loans
% R_Ds=xfs(10);        % Deposit rate
% it possible to express most of the variables of the model in close form solution (i.e. as function of the 10 variables listed above, 
% variables in 1. and parameters). 
% We can then reduce the problem to a system of 10 equations in 10 unknowns.
% Note that it is not possible to use any model equation 2 times (with each model equation we have to derive the steady state value of one single variable!)
%%%%%**********************************************************************
% Solve_SS_AdjRD.m file reduce the problem to a system of 10 equations and solves numerically for the 10 previously selected unknowns.
% initval1=[1.04 6.5 4 0.77 0.8 0.78 1.03 0.78 1.03 0.2];

%initval1=[1.07000000000000,14.8100000000000,14.5100000000000,0.620000000000000,0.250000000000000,0.870000000000000,1.03000000000000,0.870000000000000,1.03000000000000,3.97000000000000,1.00500000000000,151.060000000000,2.44000000000000,58.2000000000000,1.03000000000000,1.03000000000000,151.060000000000,1.51000000000000,-111,149.500000000000,-0.0100000000000000]
 % load initval1;
%  initval1=[1.05649,21.6402, 21.5138, 0.657728, 0.835946,0.872211,1.02623,...
%      0.872211,1.02552, 1.81102,1.00503, 129.816, 12.3594, 86.7923, 1.02623, 1.02597, 129.816 2.07706,-57.2049 127.739, -0.0390648];
%initval1=[1.05595 2.17479 1.66523  0.662915 0.773183  0.87221  1.02604 0.87221  1.02579 1.81104  1.00503  13.9905  6.07555 29.6125 1.02604 1.02579 13.9905  0.223848  -6.16514  13.7667  -0.0390162];
pu=10e-4;
load initval1.mat;
%initval1=[1.04 3.59 2.04 0.56 0.95 0.96 1.04 1.04 0.14];
options = optimset('TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',9999,'MaxIter',9999,'Display','off');
%options = optimset('TolFun',1e-6);
%[xfs fsfv,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB]=lsqnonlin('Solve_SS_AdjRD',initval1,[],[],options,psib,nu,pu,zeta_m,zeta_F,def_rate_ss,tau_H,tau_F,rp,rpe,pp,hab,delta_H,a_e,a_s,a_b,q_Ks,As,PIs,PHs,R_Hs,R_DDs,alphaa,delta_K,betta_m,betta_s,phi_F,phi_H,phi,epsilonH,epsilonF,mu_m,mu_e,mu_B,sigma_e1,sigma_m1,sigma_B,varphi_s,varphi_m,v_s,v_m,chi_b,chi_e,eta,omikronH,omikronF);
[xfs fsfv]=fsolve('Solve_SS_AdjRD',initval1,options,psib,nu,pu,zeta_m,zeta_F,def_rate_ss,tau_H,tau_F,rp,rpe,pp,hab,delta_H,a_e,a_s,a_b,q_Ks,As,PIs,PHs,R_Hs,R_DDs,alphaa,delta_K,betta_m,betta_s,phi_F,phi_H,phi,epsilonH,epsilonF,mu_m,mu_e,mu_B,sigma_e1,sigma_m1,sigma_B,varphi_s,varphi_m,v_s,v_m,chi_b,chi_e,eta,omikronH,omikronF);
%obj= @(xx) Solve_SS_AdjRD(xx,psib,nu,pu,zeta1,def_rate_ss,tau,taue,rp,rpe,pp,hab,delta_H,a_e,a_s,a_b,q_Ks,As,PIs,PHs,R_Hs,R_DDs,alphaa,delta_K,betta_m,betta_s,phi_F,phi_H,phi,epsilonH,epsilonF,mu_m,mu_e,mu_B,mu_B,sigma_e1,sigma_m1,sigma_B,sigma_B,varphi_s,varphi_m,v_s,v_m,chi_b,chi_e,eta,omikronH,omikronF);
%  [xfs,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(obj,initval1,[],[],[],[],[],[],[],options)
% % 
% if EXITFLAG==1 || EXITFLAG==3
%%
% 
% * ITEM1
% * ITEM2

% initval1=xfs;
% save initval1.mat initval1;
% % end

input=xfs;

R_Ks = xfs(1);
Ls   = xfs(2);
L_ms = xfs(3);
x_es = xfs(4);
x_ms = xfs(5);
omega_bar_Bs=xfs(6);
R_ms=xfs(7);
R_tilde_Fs=xfs(8);
phibs=xfs(9);



%%%%%**********************************************************************

%%%%%**********************************************************************
% 3. The remaining variables have close form solution in terms of the model
% parameters and the 10 variables for which we solved in 1. and 2.
%%%%%**********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Structure of the file
%   (1) Guess values for the 10 unknown variables above
%   (2) Use model relationships (FOCs and mkt clearing conditions) to compute other model variables
%   (3) Solve the system of equations consisting of the remaining 10 model relationships. 
%       Think of these equations as equilibrium 'residuals': the solver adjustes the initial guesses
%       for 10 variables above so that these 'residuals' become zero (a solution is thus found).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%-------------------------------------------------------------------------
Gamma_Bs=normcdf((log(omega_bar_Bs) - sigma_B^2/2)/sigma_B) + omega_bar_Bs*(1-normcdf((log(omega_bar_Bs)+sigma_B^2/2)/sigma_B));

G_Bs   = normcdf((log(omega_bar_Bs)-sigma_B^2/2)/sigma_B);
F_pBs   = normcdf((log(omega_bar_Bs)+sigma_B^2/2)/sigma_B);

Gamma_B_1s = (normpdf((log(omega_bar_Bs)-sigma_B^2/2)/sigma_B)-omega_bar_Bs*normpdf((log(omega_bar_Bs)+sigma_B^2/2)/sigma_B))/(sigma_B*omega_bar_Bs) + (1-normcdf((log(omega_bar_Bs)+sigma_B^2/2)/sigma_B));
def_rate_Bs =def_rate_ss+ normcdf((log(omega_bar_Bs)+sigma_B^2/2)/sigma_B)*400;
av_defs=def_rate_Bs;
 R_Ds=R_DDs/(1-pp*(av_defs));
%%-------------------------------------------------------------------------
q_Hs=1;
r_Ks = R_Ks - (1-delta_K);
YKs = r_Ks/alphaa;
Ks = ((YKs)/(Ls^(1-alphaa)))^(1/(alphaa-1));
Is = delta_K*Ks;
%%-------------------------------------------------------------------------
% ENTRERENEURS
%%-------------------------------------------------------------------------
Ys = As*Ks^(alphaa)*Ls^(1-alphaa);
ws = (1-alphaa)*Ys/Ls;
omega_bar_es = x_es/R_Ks;
G_e_1s = normpdf((log(omega_bar_es)-sigma_e1^2/2)/sigma_e1 )/(sigma_e1 *omega_bar_es);
G_es = normcdf((log(omega_bar_es)-sigma_e1^2/2)/sigma_e1 );
Gamma_es   = normcdf((log(omega_bar_es) - sigma_e1^2/2)/sigma_e1 ) + omega_bar_es*(1-normcdf((log(omega_bar_es)+sigma_e1^2/2)/sigma_e1 ));
Gamma_e_1s = (normpdf((log(omega_bar_es)-sigma_e1 ^2/2)/sigma_e1 )-omega_bar_es*normpdf((log(omega_bar_es)+sigma_e1 ^2/2)/sigma_e1 ))/(sigma_e1 *omega_bar_es) + (1-normcdf((log(omega_bar_es)+sigma_e1 ^2/2)/sigma_e1 ));
F_pes = normcdf((log(omega_bar_es)+sigma_e1^2/2)/(sigma_e1));
n_es=(1-Gamma_es)*R_Ks*q_Ks*Ks;
%xi_es=-(1-Gamma_es)*R_Ks/((1-Gamma_Bs)*(Gamma_es- mu_e*G_es) *R_Ks - rho_Fs*phi_F);


%%-------------------------------------------------------------------------
% HOUSEHOLDS
%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------
% BORROWERS
%%-------------------------------------------------------------------------
omega_bar_ms = x_ms/R_Hs;
G_m_1s = normpdf((log(omega_bar_ms)-sigma_m1^2/2)/sigma_m1)/(sigma_m1*omega_bar_ms);
Gamma_ms   = normcdf((log(omega_bar_ms)- sigma_m1^2/2)/sigma_m1) + omega_bar_ms*(1-normcdf((log(omega_bar_ms)+sigma_m1^2/2)/sigma_m1));
G_ms   = normcdf((log(omega_bar_ms)-sigma_m1^2/2)/sigma_m1);
Gamma_m_1s = (normpdf((log(omega_bar_ms)-sigma_m1^2/2)/sigma_m1)-omega_bar_ms*normpdf((log(omega_bar_ms)+sigma_m1^2/2)/sigma_m1))/(sigma_m1*omega_bar_ms) + (1-normcdf((log(omega_bar_ms)+sigma_m1^2/2)/sigma_m1));
F_pis = normcdf((log(omega_bar_ms)+sigma_m1^2/2)/(sigma_m1));


%R_ms=R_mis;
R_tilde_Hs=R_ms*((Gamma_ms - mu_m*G_ms)*R_Hs)/(x_ms);
%R_ms=1/((1-F_pis)*betta_m);
%--------------------------------
UL_m_1s=varphi_m*L_ms^(eta);
Lambda_ms=UL_m_1s/ws;
UC_m_1s=Lambda_ms;%10
C_ms=(1/UC_m_1s)/(1-hab);
%xi_ms=Lambda_ms/(rho_Hs*phi_H);

xi_ms=(Lambda_ms-((1-F_pis)*betta_m*Lambda_ms*(R_ms)))/(((R_ms)*(1-(1-F_pis)*(1-rp))));

ZZHms=betta_m*v_m/(Lambda_ms*(q_Hs) - betta_m*Lambda_ms*(1-G_ms)*R_Hs*q_Hs-xi_ms*epsilonH*q_Hs*delta_H);


%%ZZHms =betta_m*v_m/(Lambda_ms*(1-betta_m*(1-Gamma_ms)*R_Hs)-xi_ms*((1-Gamma_Bs)*(Gamma_ms- mu_m*G_ms)*R_Hs));

%%ZZHms =betta_m*v_m/(Lambda_ms*(1-betta_m*(1-Gamma_ms)*R_Hs)-xi_ms*((1-Gamma_Bs)*(Gamma_ms- mu_m*G_ms)*R_Hs)-xi_bs*(epsilon*q_Hs));
%%betta_m*UH_m_1(1) - Lambda_m*(q_H) + betta_m*Lambda_m(1)*(1-Gamma_m(1))*R_H(1)*q_H +xi_m*((1-Gamma_B(1))*(Gamma_m(1)- m_m*G_m(1))*R_H(1))*q_H-xi_m*zeta*((R_m-steady_state(R_m)))*R_m/H_m+xi_b*(epsilon*(q_H))= 0; 

%-------------------------------------------------------------------------
% SAVERS
%%-------------------------------------------------------------------------
L_ss=Ls-L_ms;
UL_s_1s=varphi_s*L_ss^(eta);%11
Lambda_ss=UL_s_1s/ws;
UC_s_1s = Lambda_ss;%9UC
C_ss=(1/UC_s_1s)/(1-hab);
Cs = C_ss + C_ms;
ZZHss =betta_s*v_s/(Lambda_ss*(1-(1-delta_H)*betta_s));

q_Hs = 1;
Hs=(ZZHss+ZZHms);
IHs = delta_H*Hs;
H_ss=ZZHss;
H_ms=ZZHms;
b_ms=x_ms*(H_ms*q_Hs)/(R_ms);
%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------



%Tr_Hs = (omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*R_tilde_Hs*(H_ms*q_Hs*x_ms)/R_ms;
%VVe=(omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*R_tilde_Fs*(q_Ks*Ks);
%W_es=((1-Gamma_es)*R_Ks*q_Ks*Ks-a_e*Tr_Hs-a_e*VVe)/(1-a_e*((omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*R_tilde_Fs*(1-chi_e)));
%n_es = (1-chi_e)*W_es;
%b_es=q_Ks*Ks-n_es;
%Tr_Fs = (omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*R_tilde_Fs*(q_Ks*Ks-(1-chi_e)*W_es);
%Trs =Tr_Fs + Tr_Hs;
W_es=((1-Gamma_es)*R_Ks*q_Ks*Ks);

n_es = (1-chi_e)*W_es;
b_es=q_Ks*Ks-n_es;
R_Fs= (x_es*(q_Ks*Ks)/(b_es));
xi_es=((R_Fs)*q_Ks*(1-F_pes)-((((1-G_es)*R_Ks*q_Ks))*betta_s))/(epsilonF*q_Ks*delta_K-((1-(1-F_pes)*(1-rpe))*R_Fs));
%Trs = (omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*((((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms))+((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es))));
%Tr_Hs = (omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*(((1-F_pis)*b_ms*(R_ms)+G_ms*(1 - mu_m)*(b_ms*(R_ms)/omega_bar_ms)));
%Tr_Fs = (omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*(R_Fs)/omega_bar_es)));
Trs=(omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*((((1-F_pis)*b_ms*R_ms+(1-F_pes)*b_es*R_Fs+...
G_ms*(1 - mu_m)*(b_ms*(R_ms)/omega_bar_ms)+G_es*(1 - mu_e)*(b_es*(R_Fs)/omega_bar_es))));

%Trs = (omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*((((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms))+((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es))));


%%-------------------------------------------------------------------------
%E_Fs = phi_F*(q_Ks*Ks - n_es);
%W_bs = (phi_F*(q_Ks*Ks - (1-chi_e)*W_es) + phi_H*b_ms)/(1-chi_b);
%R_Ds=R_DDs;
%%rhos=(((1-G_Bs)*((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms)))+(1-G_Bs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es)))-(1-F_pBs)*(b_ms*(1-phi_H))*R_Ds-(1-F_pBs)*(b_es*(1-phi_F))*R_Ds)/(b_es*phi_F+b_ms*phi_H);
% Pr_Hs =(1-G_Bs)*(((1-F_pis)*b_ms*(R_ms)+G_ms*(1 - mu_m)*(b_ms*(R_ms)/omega_bar_ms)))+(1-G_Bs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*(R_Fs)/omega_bar_es)))-(1-F_pBs)*(b_ms*(1-phibs))*R_Ds-(1-F_pBs)*(b_es*(1-phibs))*R_Ds;

% W_bs = Pr_Hs;

% C_bs=chi_b*W_bs;
%C_bs=W_bs-(phi_F*(q_Ks*Ks - (1-chi_e)*W_es) + phi_H*b_ms);
C_es=chi_e*W_es;
% Ds = (ws*L_ss - delta_H*H_ss -C_ss + C_bs + C_es - a_s*Trs+PIs+PHs)/(1-R_DDs);
% n_bs = (b_ms + (q_Ks*Ks - n_es)-Ds);
R_Fs= (x_es*(q_Ks*Ks)/(q_Ks*Ks-n_es));
omega_bar_ms= x_ms/R_Hs;
PHI_m = (1-Gamma_ms)*R_Hs*q_Hs*H_ms;



Pr_Hs = (1-G_Bs)*(((1-F_pis)*b_ms*(R_ms)+G_ms*(1 - mu_m)*(b_ms*(R_ms)/omega_bar_ms)))+(1-G_Bs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*(R_Fs)/omega_bar_es)))-(1-F_pBs)*(b_ms*(1-phibs))*R_Ds-(1-F_pBs)*(b_es*(1-phibs))*R_Ds;
W_bs= Pr_Hs;%14
C_bs	= chi_b*W_bs; %15
 Ds = (ws*L_ss - delta_H*H_ss -C_ss + C_bs + C_es - a_s*Trs+PIs+PHs)/(1-R_DDs); %16
 n_bs =(b_ms + (q_Ks*Ks - n_es)-Ds);
%A1s=((betta_s*zeta1)*(Lambda_ss*((1-F_pBs)*(R_Ds)+nu*(phibs/phi)^(1-psib)/(b_ms+b_es))*((R_ms)^(tau))*b_ms))/(1-(betta_s*zeta1));

 A1s	 = ((betta_s*zeta_m)*(Lambda_ss*((1-F_pBs)*(R_Ds)+nu*(phibs/phi_H)^(1-psib)/(b_ms+b_es))*((R_ms)^(tau_H))*b_ms))/(1-(betta_s*zeta_m));
 B1s=((betta_s*zeta_m)*(Lambda_ss*((1-G_Bs)*(((1-F_pis))+(G_ms*(1 - mu_m)/omega_bar_ms)))*((((R_ms)^(tau_H))*b_ms))))/(1-(betta_s*zeta_m));
C1s      = ((betta_s*zeta_F)*(Lambda_ss*((1-F_pBs)*(R_Ds)+nu*(phibs/phi_F)^(1-psib)/(b_ms+b_es))*((R_Fs)^(tau_F))*b_es))/(1-(betta_s*zeta_F));
D1s=(betta_s*zeta_F)*(Lambda_ss*((1-G_Bs)*(((1-F_pes))+(G_es*(1 - mu_e)/omega_bar_es)))*((R_Fs)^(tau_F)*((b_es))))/(1-(betta_s*zeta_F));
R_mis=(tau_H/(tau_H-1))*(A1s)/(B1s);
R_Fis=(tau_F/(tau_F-1))*(C1s)/(D1s);
%D1s=((betta_s*zeta)*(Lambda_ss*((1-G_Bs)*(((1-F_pes))+(G_es*(1 - mu_e)/omega_bar_es)))*((b_es)/infs1)))/(1-(betta_s*zeta));
%R_mis=R_ms;
%R_Fis=R_Fs;
%R_mis=(((A1s)/(B1s))/((C1s)/(D1s)))*R_Fis;
%R_mis=((tau/(tau-1))*(A1s)/(B1s)-(tau/(tau-1))*(C1s)/(D1s))+R_Fis;

% R_mis=(tau/(tau-1))*(A1s)/(B1s);
% R_Fis=(tau/(tau-1))*(C1s)/(D1s);

%%-betta_s*(((1-G_Bs)*(((1-F_pis)*(1-tau)+G_ms*(1 - mu_m)*(((1-tau))/omega_bar_ms))))+(1-F_pBs)*((1-phi_H))*(R_D/R_m)*tau+(phi_H*tau*rho/R_m))=xi_b*phi_H*tau/R_m=0;

%xi_bs=((1-G_Bs)*((1-F_pis*R_ms+G_ms*(1 - mu_m)*(R_ms/omega_bar_ms))));

%R_D=xi_bs/(1-F_pB(+1))*((1-phi_H));

% av_defs = ((1-phi_F)*b_es*def_rate_Bs + (1-phi_H)*b_ms*def_rate_Bs)/Ds;
def_rate_ms =def_rate_ss+ normcdf((log(omega_bar_ms)+sigma_m1^2/2)/sigma_m1,0,1)*400;
def_rate_es = def_rate_ss+normcdf((log(omega_bar_es)+sigma_e1^2/2)/sigma_e1,0,1)*400;

UC_ss = log(C_ss*(1-hab));%1
UC_ms = log(C_ms*(1-hab));%2
UL_ss = varphi_s*L_ss^(1+eta)/(1+eta);%3
UL_ms = varphi_m*L_ms^(1+eta)/(1+eta);%4
UH_ss = v_s*log(H_ss);%5
UH_ms = v_m*log(H_ms);%6
Util_ss = UC_ss - UL_ss + UH_ss;%7
Util_ms = UC_ms - UL_ms + UH_ms;%8

UH_m_1s=v_m/H_ms;
UH_s_1s=v_s/H_ss;

%%-------------------------------------------------------------------------
Welf_ss=Util_ss/(1-betta_s);
Welf_ms=Util_ms/(1-betta_m);
UT_b=n_bs;
UT_e=n_es;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  HOUSE-KEEPING STUFF BELOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
av_def=av_defs;
b_e=b_es;
A = 1; EJ=1; EK=1; EH=1;
 ESe= 1; ESm= 1; ESB= 1;
EWb=1;%ERW=1;
EWe=1;
EdH=0; EdK=0;markup_m=1;markup_F=1;EC=1;ECAB=1;EL=1;
EbH=1;EbF=1;
ELTVH=1;
ELTVF=1;
b_m = b_ms;
C_m = C_ms; C_s = C_ss; C = Cs;
D = Ds;
def_rate_e = def_rate_es;
def_rate_m = def_rate_ms;
def_rate_B = def_rate_Bs;
def_rate_B = def_rate_Bs;
%E_F = E_Fs;
G_e_1 = G_e_1s; G_e = G_es;
G_B = G_Bs;F_pB=F_pBs; G_B = G_Bs;F_pB = F_pBs;
g_I_1 = g_I_1s; g_I = g_Is;
g_H=g_Hs; g_H_1=g_H_1s;
G_m_1 = G_m_1s; G_m = G_ms;
Gamma_e = Gamma_es; Gamma_e_1 = Gamma_e_1s;
Gamma_B = Gamma_Bs; Gamma_B_1 = Gamma_B_1s;
Gamma_m = Gamma_ms; Gamma_m_1 = Gamma_m_1s;F_pi=F_pis;F_pe=F_pes;
H_m = H_ms; H_s = H_ss; H = Hs; IH=IHs;
I = Is;
K = Ks;
L_m = L_ms; L_s = L_ss; L = Ls;
Lambda_m = Lambda_ms; Lambda_s = Lambda_ss;
n_b = n_bs; n_e = n_es;
omega_bar_e = omega_bar_es;
omega_bar_B = omega_bar_Bs;
omega_bar_m = omega_bar_ms;
PI = PIs; PH=PHs;
q_H = q_Hs; q_K = q_Ks;
R_D = R_Ds; R_DD=R_DDs; R_F = R_Fs; R_H = R_Hs;
R_K = R_Ks; r_K = r_Ks;
R_m = R_ms;
R_tilde_F = R_tilde_Fs;
R_tilde_H = R_tilde_Hs;
%rho_F = rho_Fs;
%rho_H = rho_Hs;
%Tr_H = Tr_Hs; 
%Tr_F = Tr_Fs; 
Tr = Trs;
UC_m_1 = UC_m_1s;
UC_m = UC_ms;
UC_s_1 = UC_s_1s;
UC_s = UC_ss;
UH_m_1 = UH_m_1s;
UH_m = UH_ms;
UH_s_1 = UH_s_1s;
UH_s = UH_ss;
UL_m_1 = UL_m_1s;
UL_m = UL_ms;
UL_s_1 = UL_s_1s;
UL_s= UL_ss;
Util_m = Util_ms;Util_s = Util_ss;
W_b = W_bs;Pr_H=Pr_Hs; W_e = W_es;
w = ws;
x_e = x_es; x_m = x_ms;
xi_e = xi_es; xi_m = xi_ms;
Y = Ys;
A1=A1s;B1=B1s;C1=C1s;D1=D1s;R_mi=R_mis;R_Fi=R_Fis;phib=phibs;
%%-------------------------------------------------------------------------
% All _obs variables are in deviation from steady state. Hence all are zero
Y_obs = 0;
R_D_obs = 0;
R_m_obs = 0;
R_H_obs = 0;
R_F_obs = 0;
H_m_obs = 0;
H_s_obs = 0;
b_m_obs = 0;
C_obs = 0;
C_m_obs = 0;
C_s_obs = 0;
D_obs = 0;
E_F_obs = 0;
I_obs = 0;
K_obs = 0;
L_obs = 0;
L_m_obs = 0;
L_s_obs = 0;
n_b_obs = 0;
n_e_obs = 0;
q_H_obs = 0;
q_K_obs = 0;
r_K_obs = 0;
R_K_obs = 0;
R_tilde_F_obs = 0;
R_tilde_H_obs = 0;
rho_Fs_obs = 0;
rho_Hs_obs = 0;
Tr_obs = 0;
Tr_F_obs = 0;
Tr_H_obs = 0;
w_obs = 0;
x_e_obs = 0;
x_m_obs = 0;
b_e_obs = 0;

costs = mu_e*G_es*R_Ks*q_Ks*Ks ...
    + mu_m*G_ms*R_Hs*q_Hs*H_ms...
    + mu_B*G_Bs*(Gamma_es - mu_e*G_es)*R_Ks*Ks ...
    + mu_B*G_Bs*(Gamma_ms - mu_m*G_ms)*R_Hs*q_Hs*H_ms...
    + R_Ds*pp*(av_defs)*Ds;

Y_net = Cs + chi_b*W_b + chi_e*W_e + Is + IHs;
Y_net_obs = 0;
Y_net_2 = Ys - costs;
Y_net_2_obs = 0;
GDP = Ys + UL_s_1*H_s/Lambda_s + UL_m_1*H_m/Lambda_m;
H = H_s + H_m;
H_obs = 0;
IH_obs = 0;
GDP_obs = 0;
res_H = 0;
m_e = mu_e;
m_m = mu_m;
R_D_rf = R_D;
RDsp =400*(R_D-R_DD);
W_m = PHI_m + w*L_m;
W_m_obs = 0;
b_tot=b_e+b_m;
b_tot_obs=0;
res_chk = (R_D*pp*av_def*D+C + chi_b*W_b + chi_e*W_e + I + IH + (m_e*G_e*R_K*q_K*K + m_m*G_m*R_H*q_H*H_m+ mu_B*G_B*(Gamma_e - m_e*G_e)*R_K*q_K*K + mu_B*G_B*(Gamma_m - m_m*G_m)*R_H*q_H*H_m))/Y;



%=====
%rho_Fs=(1-Gamma_Bs)*R_tilde_Fs/phi_Fs;
%rho_Hs=(1-Gamma_Bs)*R_tilde_Hs/phi_Hs;
roe_Bs=((1-Gamma_Bs)*(b_ms*R_tilde_Hs+b_es*R_tilde_Fs)/(b_ms+b_es))/phibs;




DY_NET=0;
welfare_=(C_ms*Welf_ms+C_ss*Welf_ss)/(C_ms+C_ss);
welfare_obs=0;
DY_NET_OBS=0;
labor_growth_data=0.095;
credit_gap_obs=0;
bsp_H_obs=0;
bsp_F_obs=0;
%===MEASUREMENT EQUATION VARIABLES
dy_data=gamma_y;
dq_H_data=gamma_dq_H;
dbe_data=gamma_dbe;
dbm_data=gamma_dbm;
db_total_data=5.73;
bank_rate_data=100*(R_D-1);
b_to_Y_data=0;
d_b_to_Y_data=0;
int_rate_HH_data=100*(R_m-1);
int_rate_business_data=100*(R_F-1);
pb_ratio_data=pb_mean;

bsp_H_data=0;
bsp_F_data=0;
% bsp_H_data=400*(R_tilde_Hs-R_Ds)+gamma_bspH;
%bsp_F_data=400*(R_tilde_Fs-R_Ds);
%bsp_H_data=400*(R_tilde_F-R_D);
%bsp_H_data=gamma_bspH;
% bsp_F_data=gamma_bspH;
dw_data=gamma_w;
dinve_data=gamma_inve;
dc_data=gamma_c;
CR_data=100*phi_Hs;
%HH_debt_to_GDP=0;
%HH_debt_to_GDP=(b_ms/Y_net);
%private_debt_to_GDP=(b_es/Y_net);
%roes=rho_Fs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ENTIRE SS TO BE LOADED INTO DYNARE    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ys=[b_e%1
        g_H%2
        g_H_1%3
        IH%4
        PH%5
        A%6
        b_m%7
        C%8
        C_m%9
        C_s%10
        D%11
        def_rate_e%12
        def_rate_m%13
        def_rate_B%14
        %def_rate_B%15
        epsilonH%16
        epsilonF%17
        G_e%18
        G_e_1%19
        %G_B%20
        F_pB%21
        G_B%22
        %F_pB%23
        g_I%24
        g_I_1%25
        G_m%26
        G_m_1%27
        Gamma_e%28
        Gamma_e_1%29
        %Gamma_B%30
        %Gamma_B_1%31
        Gamma_B%32
        Gamma_B_1%33
        Gamma_m%34
        Gamma_m_1%35
        F_pi%36
        F_pe%37
        H%38
        H_m%39
        H_s%40
        I%41
        K%42
        L%43
        L_m%44
        L_s%45
        Lambda_m%46
        Lambda_s%47
        n_b%48
        n_e%49
        omega_bar_e%50
        %omega_bar_B%51
        omega_bar_B%52
        omega_bar_m%53
        PI%54
        q_H%55
        q_K%56
        R_D%57
        R_F%58
        R_H%59
        r_K%60
        R_K%61
        R_m%62
        R_tilde_F%63
        R_tilde_H%64
        %Tr_H%65
        %Tr_F%66
        Tr%67
        UC_m%68
        UC_m_1%69
        UC_s%70
        UC_s_1%71
        UH_m%72
        UH_m_1%73
        UH_s%74
        UH_s_1%75
        UL_m%76
        UL_m_1%77
        UL_s%78
        UL_s_1%79
        Util_m%80
        Util_s%81
        w%82
        W_b%83
        Pr_H%84
        W_e%85
        x_e%86
        x_m%87
        xi_e%88
        xi_m%89
        Y%90
        Y_obs%91
        R_D_obs%92
        R_m_obs%93
        R_H_obs%94
        R_F_obs%95
        H_m_obs%96
        H_s_obs%97
        b_m_obs%98
        C_obs%99
        C_m_obs%100
        C_s_obs%101
        D_obs%102
        I_obs%103
        K_obs%104
        L_obs%105
        L_m_obs%106
        L_s_obs%107
        n_b_obs%108
        n_e_obs%109
        q_H_obs%110
        q_K_obs%111
        r_K_obs%112
        R_K_obs%113
       % R_tilde_F_obs%114
       % R_tilde_H_obs%115
       % Tr_obs%116
        w_obs%117
        x_e_obs%118
        x_m_obs%119
        Welf_ss%120
        Welf_ms%121
        UT_e%122
        UT_b%123
        EJ%124
        EK%125
        EH%126
        ESe%127
        ESm%128
        %ESH%129
        %ESF%130
        ESB
        EWe%131
        EWb%132
        EdH%133
        EdK%134
        markup_m
        markup_F
        EC
        ECAB
        EL
        EbH
        EbF
        ELTVH
        ELTVF
        phi_F%135
        phi_H%136
        phi%137
        b_e_obs%138
      400*  (R_tilde_F-R_D)%139
      400*  (R_tilde_H-R_D)%140
        Y_net%141
        Y_net_obs%142
        Y_net_2%143
        Y_net_2_obs%144
        res_chk%145
        delta_H%146
        delta_K%147
        GDP%148
        GDP_obs%149
        H_obs%150
        IH_obs%151
        res_H%152
        m_e%153
        m_m%154
        av_def%155
        RDsp%156
        R_DD%157
        W_m%158
        W_m_obs%159
        b_tot%160
        b_tot_obs%161
        C_bs%162
        C_es%163
        A1%164
        B1%165
        C1%166
        D1%167
        R_mi%168
        R_Fi%169
        phib%170
        %rho_Fs
         %rho_Hs
         %rho_Fs_obs
         %rho_Hs_obs
         roe_Bs
         DY_NET
         DY_NET_OBS
         

         
         credit_gap_obs
         bsp_H_obs
         bsp_F_obs
        dy_data
        dq_H_data
        dbe_data
        dbm_data
        db_total_data
        bank_rate_data
        b_to_Y_data
        d_b_to_Y_data
        int_rate_HH_data
        int_rate_business_data
       %bsp_H_data
       % bsp_F_data
        dw_data
        dinve_data
        dc_data
        CR_data
        welfare_
        welfare_obs
        labor_growth_data
        pb_ratio_data

    ];
