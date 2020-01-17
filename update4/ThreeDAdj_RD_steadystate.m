%%%%%**********************************************************************
% Codes to solve the 3D MODEL
% UPON USE PLEASE CITE: "Capital Regulation in a Macroeconomic Model with Three Layers of Default"
% by  Laurent Clerc, Alexis Derviz, Caterina Mendicino, Stephane Moyen, Kalin Nikolov, 
% Livio Stracca, Javier Suarez, Alexandros P. Vardoulakis, 
% International Journal of Central Banking, June 2015, Pages 9-63.
%%%%%**********************************************************************
% THIS PROGRAM WAS TESTED WITH (AND REQUIRES KNOWLEDGE OF)
% 1) DYNARE 4.4.3 (http://www.cepremap.cnrs.fr/dynare/)
% 2) MATLAB R2012b (http://www.mathworks.com/)
%%%%%**********************************************************************
% Contact persons:
% CATERINA MENDICINO: caterina.mendicino1@ecb.int
% KALIN NIKOLOV     : kalin.Nikolov@ecb.int
% DOMINIK SUPERA    : dominik.supera@gmail.com
%%%%%**********************************************************************
% Latest update: 22/02/2016 by DOMINIK SUPERA
%%%%%**********************************************************************

%%%%%**********************************************************************
% This file solves for the steady state of the model (deterministic equilibrium in absence of shocks).
% The steady state file contains:
% 1. close form solution for a first set of endogenous variables (easily
% derived from a subset of model's equations)
% 2. numerical solution for a subset of 10 variables:
% Solve_SS_AdjRD.m is a file that reduces the problem to a system of 10 equations and solves numerically for the 10 previously selected unknow
% 3. derivations of the remaining variables in terms of the model parameters and the variables for which we solved in 1. and 2.
%%%%%**********************************************************************

function [ys,check,input] = ThreeDAdj_RD_steadystate(junk,ys)
format long
check=0;

%%-------------------------------------------------------------------------
%Load mat file with parameter values structure
load par_ThreeDAdj_RD_mod;
filename = 'par_ThreeDAdj_RD_mod.mat';
m = matfile(filename,'Writable',true);

m.pp=0;
m.phis=0.2;
m.psib=10;
m.nu=0.3;
m.chi_b=0.1;
m.chi_e=0.1;
m.alphaa=0.2;
m.pp=0;
m.pu=0.0001;
m.rp=0.035;
m.rpe=0.04;
m.delta_H=0.035;
m.delta_K=0.035;
m.tau=50;
m.taue=50;
m.phi_Fs = 0.2;
m.phi_Hs = 0.2;
m.zeta1=0.02;
m.betta_m=0.97;
m.epsilonH1s=0.8;
m.epsilonF1s=0.8;
m.v_m=0.25;
m.v_s=0.25;
m.mu_H=0.3;
m.mu_m=0.25;
m.mu_e=0.25;
m.sigma_m1=0.12;
m.sigma_e1=0.2;
m.sigma_epsiRW=25;
m.sigma_H=0.15;
m.sigma_F=0.15;

%%-------------------------------------------------------------------------

%%%%%**********************************************************************   
% 1) A first set of endogenous variables has close form solution and 
% can be easily derived from a subset of equations
%%%%%**********************************************************************
%ERW=1;
phi_F=phi_Fs;
phi_H=phi_Hs;
phi=phis;
infs1=1;
%ELTV=1;
epsilonH=epsilonH1s;
epsilonF=epsilonF1s;
R_DDs = (1/betta_s)*infs1;
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
% omega_bar_Hs=xfs(6); % Mortgage bank default cut off 
% R_tilde_Hs=xfs(7);   % Return on portfolio of mortgage loans
% omega_bar_Fs=xfs(8); % Corporate bank default cut off
% R_tilde_Fs=xfs(9);   % Return on portfolio of corporate loans
% R_Ds=xfs(10);        % Deposit rate
% it possible to express most of the variables of the model in close form solution (i.e. as function of the 10 variables listed above, 
% variables in 1. and parameters). 
% We can then reduce the problem to a system of 10 equations in 10 unknowns.
% Note that it is not possible to use any model equation 2 times (with each model equation we have to derive the steady state value of one single variable!)
%%%%%**********************************************************************
% Solve_SS_AdjRD.m file reduce the problem to a system of 10 equations and solves numerically for the 10 previously selected unknowns.

options = optimset('TolFun',1e-28,'MaxFunEvals',3000);
[xfs fsfv]=fsolve('Solve_SS_AdjRD',[input_1,input_2,input_3,input_4,input_5,input_6,input_7,input_8,input_9,input_10],options,psib,nu,pu,zeta1,tau,taue,rp,rpe,pp,hab,delta_H,a_e,a_s,a_b,q_Ks,As,PIs,PHs,R_Hs,R_DDs,alphaa,delta_K,betta_m,betta_s,phi_F,phi_H,phi,epsilonH,epsilonF,mu_m,mu_e,mu_F,mu_H,sigma_e1,sigma_m1,sigma_F,sigma_H,varphi_s,varphi_m,v_s,v_m,chi_b,chi_e,eta);
input_1 = xfs(1);
input_2   = xfs(2);
input_3 = xfs(3);
input_4 = xfs(4);
input_5 = xfs(5);
input_6 =xfs(6);
input_7 =xfs(7);
input_8 =xfs(8);
input_9 =xfs(9);
input_10 =xfs(10);

input=xfs;

R_Ks = xfs(1);
Ls   = xfs(2);
L_ms = xfs(3);
x_es = xfs(4);
x_ms = xfs(5);
omega_bar_Hs=xfs(6);
R_ms=xfs(7);
omega_bar_Fs=xfs(8);
R_tilde_Fs=xfs(9);
phibs=xfs(10);

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
Gamma_Hs=normcdf((log(omega_bar_Hs)-sigma_H^2/2)/sigma_H) + omega_bar_Hs*(1-normcdf((log(omega_bar_Hs)+sigma_H^2/2)/sigma_H));
Gamma_Fs=normcdf((log(omega_bar_Fs) - sigma_F^2/2)/sigma_F) + omega_bar_Fs*(1-normcdf((log(omega_bar_Fs)+sigma_F^2/2)/sigma_F));
G_Hs   = normcdf((log(omega_bar_Hs)-sigma_H^2/2)/sigma_H);
F_pHs   = normcdf((log(omega_bar_Hs)+sigma_H^2/2)/sigma_H);

Gamma_H_1s = (normpdf((log(omega_bar_Hs)-sigma_H^2/2)/sigma_H)-omega_bar_Hs*normpdf((log(omega_bar_Hs)+sigma_H^2/2)/sigma_H))/(sigma_H*omega_bar_Hs) + (1-normcdf((log(omega_bar_Hs)+sigma_H^2/2)/sigma_H));
def_rate_Hs = normcdf((log(omega_bar_Hs)+sigma_H^2/2)/sigma_H)*400;
G_Fs   = normcdf((log(omega_bar_Fs)-sigma_F^2/2)/sigma_F);
F_pFs   = normcdf((log(omega_bar_Fs)+sigma_F^2/2)/sigma_F);

Gamma_F_1s = (normpdf((log(omega_bar_Fs)-sigma_F^2/2)/sigma_F)-omega_bar_Fs*normpdf((log(omega_bar_Fs)+sigma_F^2/2)/sigma_F))/(sigma_F*omega_bar_Fs) + (1-normcdf((log(omega_bar_Fs)+sigma_F^2/2)/sigma_F));
def_rate_Fs = normcdf((log(omega_bar_Fs)+sigma_F^2/2)/sigma_F)*400;
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
%xi_es=-(1-Gamma_es)*R_Ks/((1-Gamma_Fs)*(Gamma_es- mu_e*G_es) *R_Ks - rho_Fs*phi_F);


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
R_tilde_Hs=R_ms*((Gamma_ms - mu_m*G_ms)*R_Hs)*infs1/(x_ms);
%R_ms=1/((1-F_pis)*betta_m);
%--------------------------------
UL_m_1s=varphi_m*L_ms^(eta);
Lambda_ms=UL_m_1s/ws;
UC_m_1s=Lambda_ms;
C_ms=(1/UC_m_1s)/(1-hab);
%xi_ms=Lambda_ms/(rho_Hs*phi_H);

xi_ms=(Lambda_ms-((1-F_pis)*betta_m*Lambda_ms*(R_ms/infs1)))/(((R_ms/infs1)*(1-(1-F_pis)*(1-rp))));

ZZHms=betta_m*v_m/(Lambda_ms*(q_Hs) - betta_m*Lambda_ms*(1-G_ms)*R_Hs*q_Hs-xi_ms*epsilonH*q_Hs*delta_H);


%%ZZHms =betta_m*v_m/(Lambda_ms*(1-betta_m*(1-Gamma_ms)*R_Hs)-xi_ms*((1-Gamma_Hs)*(Gamma_ms- mu_m*G_ms)*R_Hs));

%%ZZHms =betta_m*v_m/(Lambda_ms*(1-betta_m*(1-Gamma_ms)*R_Hs)-xi_ms*((1-Gamma_Hs)*(Gamma_ms- mu_m*G_ms)*R_Hs)-xi_bs*(epsilon*q_Hs));
%%betta_m*UH_m_1(1) - Lambda_m*(q_H) + betta_m*Lambda_m(1)*(1-Gamma_m(1))*R_H(1)*q_H +xi_m*((1-Gamma_H(1))*(Gamma_m(1)- m_m*G_m(1))*R_H(1))*q_H-xi_m*zeta*((R_m-steady_state(R_m)))*R_m/H_m+xi_b*(epsilon*(q_H))= 0; 

%-------------------------------------------------------------------------
% SAVERS
%%-------------------------------------------------------------------------
L_ss=Ls-L_ms;
UL_s_1s=varphi_s*L_ss^(eta);
Lambda_ss=UL_s_1s/ws;
UC_s_1s = Lambda_ss;
C_ss=(1/UC_s_1s)/(1-hab);
Cs = C_ss + C_ms;
ZZHss =betta_s*v_s/(Lambda_ss*(1-(1-delta_H)*betta_s));

q_Hs = 1;
Hs=(ZZHss+ZZHms);
IHs = delta_H*Hs;
H_ss=ZZHss;
H_ms=ZZHms;
b_ms=x_ms*(H_ms*q_Hs)/(R_ms/infs1);
%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------



%Tr_Hs = (omega_bar_Hs - Gamma_Hs + mu_H*G_Hs)*R_tilde_Hs*(H_ms*q_Hs*x_ms)/R_ms;
%VVe=(omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*R_tilde_Fs*(q_Ks*Ks);
%W_es=((1-Gamma_es)*R_Ks*q_Ks*Ks-a_e*Tr_Hs-a_e*VVe)/(1-a_e*((omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*R_tilde_Fs*(1-chi_e)));
%n_es = (1-chi_e)*W_es;
%b_es=q_Ks*Ks-n_es;
%Tr_Fs = (omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*R_tilde_Fs*(q_Ks*Ks-(1-chi_e)*W_es);
%Trs =Tr_Fs + Tr_Hs;
W_es=((1-Gamma_es)*R_Ks*q_Ks*Ks);

n_es = (1-chi_e)*W_es;
b_es=q_Ks*Ks-n_es;
R_Fs= (x_es*(q_Ks*Ks)/(b_es))*infs1;
xi_es=((R_Fs/infs1)*q_Ks*(1-F_pes)-((((1-G_es)*R_Ks*q_Ks))*betta_s))/(epsilonF*q_Ks*delta_K-((1-(1-F_pes)*(1-rpe))*R_Fs/infs1));
%Trs = (omega_bar_Hs - Gamma_Hs + mu_H*G_Hs)*((((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms))+((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es))));
Tr_Hs = (omega_bar_Hs - Gamma_Hs + mu_H*G_Hs)*(((1-F_pis)*b_ms*(R_ms/infs1)+G_ms*(1 - mu_m)*(b_ms*(R_ms/infs1)/omega_bar_ms)));
Tr_Fs = (omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*(((1-F_pes)*b_es*R_Fs/infs1+G_es*(1 - mu_e)*(b_es*(R_Fs/infs1)/omega_bar_es)));
Trs = Tr_Fs + Tr_Hs;
%Trs = (omega_bar_Hs - Gamma_Hs + mu_H*G_Hs)*((((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms))+((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es))));


%%-------------------------------------------------------------------------
%E_Fs = phi_F*(q_Ks*Ks - n_es);
%W_bs = (phi_F*(q_Ks*Ks - (1-chi_e)*W_es) + phi_H*b_ms)/(1-chi_b);
R_Ds=R_DDs;
%%rhos=(((1-G_Hs)*((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms)))+(1-G_Fs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es)))-(1-F_pHs)*(b_ms*(1-phi_H))*R_Ds-(1-F_pFs)*(b_es*(1-phi_F))*R_Ds)/(b_es*phi_F+b_ms*phi_H);
Pr_Hs =(1-G_Hs)*(((1-F_pis)*b_ms*(R_ms/infs1)+G_ms*(1 - mu_m)*(b_ms*(R_ms/infs1)/omega_bar_ms)))+(1-G_Fs)*(((1-F_pes)*b_es*R_Fs/infs1+G_es*(1 - mu_e)*(b_es*(R_Fs/infs1)/omega_bar_es)))-(1-F_pHs)*(b_ms*(1-phibs))*R_Ds-(1-F_pFs)*(b_es*(1-phibs))*R_Ds;

W_bs = Pr_Hs;

C_bs=chi_b*W_bs;
%C_bs=W_bs-(phi_F*(q_Ks*Ks - (1-chi_e)*W_es) + phi_H*b_ms);
C_es=chi_e*W_es;
Ds = (ws*L_ss - delta_H*H_ss -C_ss + C_bs + C_es - a_s*Trs+PIs+PHs)/(1-R_DDs);
n_bs = (b_ms + (q_Ks*Ks - n_es)-Ds);
R_Fs= (x_es*(q_Ks*Ks)/(q_Ks*Ks-n_es))*infs1;
omega_bar_ms= x_ms/R_Hs;
PHI_m = (1-Gamma_ms)*R_Hs*q_Hs*H_ms;

%rhos=xi_bs;
%rhos=(((1-G_Hs)*((1-F_pis)*b_ms*(R_ms/infs1)+G_ms*(1 - mu_m)*(b_ms*(R_ms/infs1)/omega_bar_ms)))+(1-G_Fs)*(((1-F_pes)*b_es*R_Fs/infs1+G_es*(1 - mu_e)*(b_es*(R_Fs/infs1)/omega_bar_es)))-(1-F_pHs)*(b_ms*(1-phi_H))*R_Ds-(1-F_pFs)*(b_es*(1-phi_F))*R_Ds)/(b_es*phi_F+b_ms*phi_H);
%xi_bs=-betta_s*(((((1-G_Hs)*(((1-F_pis)*(1-tau)+G_ms*(1 - mu_m)*((1-tau)/omega_bar_ms))))))+((R_Ds/(R_ms/infs1))*(tau)*((1-F_pHs)*(1-phi_H))+(phi_H*tau*rhos/(R_ms/infs1))*0))*(R_ms/infs1)/(phi_H*tau);


%xi_bs=betta_s*(((1-G_Hs)*(((1-F_pis)*(tau-1)*(R_ms/infs1)+G_ms*(1 - mu_m)*(tau-1)*((R_ms/infs1)/omega_bar_ms)))-(1-F_pHs)*((tau)*(1-phi_H))*R_Ds))/(phi_H*tau);
%phibs=n_bs/(b_ms+b_es);


%A1s=((betta_s*zeta)*(Lambda_ss*((1-F_pHs)*((1-phi_H))*(R_Ds/infs1)+rhos*phi_H)*((R_ms)^(tau))*b_ms))/(1-(betta_s*zeta));


A1s=((betta_s*zeta1)*(Lambda_ss*((1-F_pHs)*(R_Ds/infs1)+nu*(phibs/phi)^(1-psib)/(b_ms+b_es))*((R_ms)^(tau))*b_ms))/(1-(betta_s*zeta1));

%A1s=((betta_s*zeta)*(Lambda_ss*(((rhos*(phi_F)+(1-F_pFs)*(1-phi_F)*R_Ds)))*b_ms))/(1-(betta_s*zeta));

B1s=((betta_s*zeta1)*(Lambda_ss*((1-G_Hs)*(((1-F_pis))+(G_ms*(1 - mu_m)/omega_bar_ms)))*((((R_ms)^(tau))*b_ms)/infs1)))/(1-(betta_s*zeta1));
%B1s=((betta_s*zeta)*(Lambda_ss*((1-G_Hs)*(((1-F_pis))+(G_ms*(1 - mu_m)/omega_bar_ms)))*((b_ms)/infs1)))/(1-(betta_s*zeta));


%C1s=((betta_s*zeta)*(Lambda_ss*((1-F_pFs)*((1-phi_F))*(R_Ds/infs1)+rhos*phi_F)*((R_Fs)^(tau))*b_es))/(1-(betta_s*zeta));
C1s=((betta_s*zeta1)*(Lambda_ss*((1-F_pFs)*(R_Ds/infs1)+nu*(phibs/phi)^(1-psib)/(b_ms+b_es))*((R_Fs)^(tau))*b_es))/(1-(betta_s*zeta1));

%C1s=((betta_s*zeta)*(Lambda_ss*((1-F_pFs)*((1-phi_F))*(R_Ds/infs1)+rhos*phi_F)*b_es))/(1-(betta_s*zeta));
D1s=(betta_s*zeta1)*(Lambda_ss*((1-G_Fs)*(((1-F_pes))+(G_es*(1 - mu_e)/omega_bar_es)))*((R_Fs)^(tau)*((b_es)/infs1)))/(1-(betta_s*zeta1));

%D1s=((betta_s*zeta)*(Lambda_ss*((1-G_Fs)*(((1-F_pes))+(G_es*(1 - mu_e)/omega_bar_es)))*((b_es)/infs1)))/(1-(betta_s*zeta));
%R_mis=R_ms;
%R_Fis=R_Fs;
%R_mis=(((A1s)/(B1s))/((C1s)/(D1s)))*R_Fis;
%R_mis=((tau/(tau-1))*(A1s)/(B1s)-(tau/(tau-1))*(C1s)/(D1s))+R_Fis;

R_mis=(tau/(tau-1))*(A1s)/(B1s);
R_Fis=(tau/(tau-1))*(C1s)/(D1s);

%%-betta_s*(((1-G_Hs)*(((1-F_pis)*(1-tau)+G_ms*(1 - mu_m)*(((1-tau))/omega_bar_ms))))+(1-F_pHs)*((1-phi_H))*(R_D/R_m)*tau+(phi_H*tau*rho/R_m))=xi_b*phi_H*tau/R_m=0;

%xi_bs=((1-G_Hs)*((1-F_pis*R_ms+G_ms*(1 - mu_m)*(R_ms/omega_bar_ms))));

%R_D=xi_bs/(1-F_pH(+1))*((1-phi_H));

av_defs = ((1-phi_F)*b_es*def_rate_Fs + (1-phi_H)*b_ms*def_rate_Hs)/Ds;
def_rate_ms = normcdf((log(omega_bar_ms)+sigma_m1^2/2)/sigma_m1)*400;
def_rate_es = normcdf((log(omega_bar_es)+sigma_e1^2/2)/sigma_e1)*400;

UC_ss = log(C_ss*(1-hab));
UC_ms = log(C_ms*(1-hab));
UL_ss = varphi_s*L_ss^(1+eta)/(1+eta);
UL_ms = varphi_m*L_ms^(1+eta)/(1+eta);
UH_ss = v_s*log(H_ss);
UH_ms = v_m*log(H_ms);
Util_ss = UC_ss - UL_ss + UH_ss;
Util_ms = UC_ms - UL_ms + UH_ms;

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
ESe= 1; ESm= 1; ESH= 1; ESF= 1;
EWb=1;%ERW=1;
EWe=1;
EdH=0; EdK=0;
b_m = b_ms;
C_m = C_ms; C_s = C_ss; C = Cs;
D = Ds;
def_rate_e = def_rate_es;
def_rate_m = def_rate_ms;
def_rate_F = def_rate_Fs;
def_rate_H = def_rate_Hs;
%E_F = E_Fs;
G_e_1 = G_e_1s; G_e = G_es;
G_F = G_Fs;F_pF=F_pFs; G_H = G_Hs;F_pH = F_pHs;
g_I_1 = g_I_1s; g_I = g_Is;
g_H=g_Hs; g_H_1=g_H_1s;
G_m_1 = G_m_1s; G_m = G_ms;
Gamma_e = Gamma_es; Gamma_e_1 = Gamma_e_1s;
Gamma_F = Gamma_Fs; Gamma_F_1 = Gamma_F_1s;
Gamma_H = Gamma_Hs; Gamma_H_1 = Gamma_H_1s;
Gamma_m = Gamma_ms; Gamma_m_1 = Gamma_m_1s;F_pi=F_pis;F_pe=F_pes;
H_m = H_ms; H_s = H_ss; H = Hs; IH=IHs;
I = Is;
K = Ks;
L_m = L_ms; L_s = L_ss; L = Ls;
Lambda_m = Lambda_ms; Lambda_s = Lambda_ss;
n_b = n_bs; n_e = n_es;
omega_bar_e = omega_bar_es;
omega_bar_F = omega_bar_Fs;
omega_bar_H = omega_bar_Hs;
omega_bar_m = omega_bar_ms;
PI = PIs; PH=PHs;
q_H = q_Hs; q_K = q_Ks;
R_D = R_Ds; R_DD=R_DDs; R_F = R_Fs; R_H = R_Hs;
R_K = R_Ks; r_K = r_Ks;
R_m = R_ms;
inf1=infs1;
R_tilde_F = R_tilde_Fs;
R_tilde_H = R_tilde_Hs;
%rho_F = rho_Fs;
%rho_H = rho_Hs;
Tr_H = Tr_Hs; 
Tr_F = Tr_Fs; 
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
%rho_F_obs = 0;
%rho_H_obs = 0;
Tr_obs = 0;
Tr_F_obs = 0;
Tr_H_obs = 0;
w_obs = 0;
x_e_obs = 0;
x_m_obs = 0;
b_e_obs = 0;

costs = mu_e*G_es*R_Ks*q_Ks*Ks ...
    + mu_m*G_ms*R_Hs*q_Hs*H_ms...
    + mu_F*G_Fs*(Gamma_es - mu_e*G_es)*R_Ks*Ks ...
    + mu_H*G_Hs*(Gamma_ms - mu_m*G_ms)*R_Hs*q_Hs*H_ms...
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
res_chk = (R_D*pp*av_def*D+C + chi_b*W_b + chi_e*W_e + I + IH + (m_e*G_e*R_K*q_K*K + m_m*G_m*R_H*q_H*H_m+ mu_F*G_F*(Gamma_e - m_e*G_e)*R_K*q_K*K + mu_H*G_H*(Gamma_m - m_m*G_m)*R_H*q_H*H_m))/Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ENTIRE SS TO BE LOADED INTO DYNARE    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ys=[b_e
        g_H
        g_H_1
        IH
        PH
        A
        b_m
        C
        C_m
        C_s
        D
        def_rate_e
        def_rate_m
        def_rate_F
        def_rate_H
        epsilonH
        epsilonF
        G_e
        G_e_1
        G_F
        F_pF
        G_H
        F_pH
        g_I
        g_I_1
        G_m
        G_m_1
        Gamma_e
        Gamma_e_1
        Gamma_F
        Gamma_F_1
        Gamma_H
        Gamma_H_1
        Gamma_m
        Gamma_m_1
        F_pi
        F_pe
        H
        H_m
        H_s
        I
        K
        L
        L_m
        L_s
        Lambda_m
        Lambda_s
        n_b
        n_e
        omega_bar_e
        omega_bar_F
        omega_bar_H
        omega_bar_m
        PI
        q_H
        q_K
        R_D
        R_F
        R_H
        r_K
        R_K
        R_m
        inf1
        R_tilde_F
        R_tilde_H
        Tr_H
        Tr_F
        Tr
        UC_m
        UC_m_1
        UC_s
        UC_s_1
        UH_m
        UH_m_1
        UH_s
        UH_s_1
        UL_m
        UL_m_1
        UL_s
        UL_s_1
        Util_m
        Util_s
        w
        W_b
        Pr_H
        W_e
        x_e
        x_m
        xi_e
        xi_m
        Y
        Y_obs
        R_D_obs
        R_m_obs
        R_H_obs
        R_F_obs
        H_m_obs
        H_s_obs
        b_m_obs
        C_obs
        C_m_obs
        C_s_obs
        D_obs
        I_obs
        K_obs
        L_obs
        L_m_obs
        L_s_obs
        n_b_obs
        n_e_obs
        q_H_obs
        q_K_obs
        r_K_obs
        R_K_obs
        R_tilde_F_obs
        R_tilde_H_obs
        Tr_obs
        w_obs
        x_e_obs
        x_m_obs
        Welf_ss
        Welf_ms
        UT_e
        UT_b
        EJ
        EK
        EH
        ESe
        ESm
        ESH
        ESF
        EWe
        EWb
        EdH
        EdK
        phi_F
        phi_H
        phi
        b_e_obs
        400*(R_tilde_F-R_D)
        400*(R_tilde_H-R_D)
        Y_net
        Y_net_obs
        Y_net_2
        Y_net_2_obs
        res_chk
        delta_H
        delta_K
        GDP
        GDP_obs
        H_obs
        IH_obs
        res_H
        m_e
        m_m
        av_def
        RDsp
        R_DD
        W_m
        W_m_obs
        b_tot
        b_tot_obs
        C_bs
        C_es
        A1
        B1
        C1
        D1
        R_mi
        R_Fi
        phib
        ];
