clear;clc;close all;
syms psib nu pu zeta1 tau taue rp rpe pp hab delta_H a_e a_s a_b q_Ks As PIs PHs R_Hs R_DDs alphaa delta_K betta_m betta_s phi_F phi_H phi epsilonH epsilonF mu_m mu_e mu_F mu_H sigma_e1 sigma_m1 sigma_F sigma_H varphi_s varphi_m v_s v_m chi_b chi_e eta;
syms R_Ks Ls L_ms x_es x_ms omega_bar_Hs R_ms omega_bar_Fs R_tilde_Fs phibs R_Ds 
variables= [ R_Ks Ls L_ms x_es x_ms omega_bar_Hs R_ms omega_bar_Fs R_tilde_Fs phibs R_Ds ];

% R_Ks = x(1);
% Ls   = x(2);
% L_ms = x(3);
% x_es = x(4);
% x_ms = x(5);           		 
% omega_bar_Hs=x(6);
% R_ms=x(7);
% omega_bar_Fs=x(8);
% R_tilde_Fs=x(9);
% phibs=x(10);
% R_Ds=x(11);
R_Ds=R_DDs;
%%rhos = 1/(1-chi_b);
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


%%-----------------------------------------------------------------
r_Ks = R_Ks - (1-delta_K);
q_Hs=1;
%rhos = 1/(1-chi_b);
YKs = r_Ks/alphaa;

Ks = ((YKs)/(Ls^(1-alphaa)))^(1/(alphaa-1)); %%% CHANGED Ks = ((YKs)/(Ls^(1-alphaa)))^(alphaa-1)

Is = delta_K*Ks;

%%-------------------------------------------------------------------------
%% ENTREPRENEURS
%%-------------------------------------------------------------------------
Ys = As*Ks^(alphaa)*Ls^(1-alphaa);

ws = (1-alphaa)*Ys/Ls;

omega_bar_es = x_es/R_Ks;

G_e_1s = normpdf((log(omega_bar_es)-sigma_e1^2/2)/sigma_e1 )/(sigma_e1 *omega_bar_es);

G_es = normcdf((log(omega_bar_es)-sigma_e1^2/2)/sigma_e1 );

Gamma_es   = normcdf((log(omega_bar_es) - sigma_e1^2/2)/sigma_e1 ) + omega_bar_es*(1-normcdf((log(omega_bar_es)+sigma_e1^2/2)/sigma_e1 ));

Gamma_e_1s = (normpdf((log(omega_bar_es)-sigma_e1 ^2/2)/sigma_e1 )-omega_bar_es*normpdf((log(omega_bar_es)+sigma_e1 ^2/2)/sigma_e1 ))/(sigma_e1 *omega_bar_es) + (1-normcdf((log(omega_bar_es)+sigma_e1 ^2/2)/sigma_e1 ));


%xi_es=-(1-Gamma_es)*R_Ks/((1-Gamma_Fs)*(Gamma_es- mu_e*G_es) *R_Ks - rho_Fs*phi_F);
%%-------------------------------------------------------------------------
%% HOUSEHOLDS
%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------
%% BORROWERS
%%-------------------------------------------------------------------------
omega_bar_ms = x_ms/R_Hs;

G_m_1s = normpdf((log(omega_bar_ms)-sigma_m1^2/2)/sigma_m1)/(sigma_m1*omega_bar_ms);

Gamma_ms   = normcdf((log(omega_bar_ms)- sigma_m1^2/2)/sigma_m1) + omega_bar_ms*(1-normcdf((log(omega_bar_ms)+sigma_m1^2/2)/sigma_m1));

G_ms   = normcdf((log(omega_bar_ms)-sigma_m1^2/2)/sigma_m1);

Gamma_m_1s = (normpdf((log(omega_bar_ms)-sigma_m1^2/2)/sigma_m1)-omega_bar_ms*normpdf((log(omega_bar_ms)+sigma_m1^2/2)/sigma_m1))/(sigma_m1*omega_bar_ms) + (1-normcdf((log(omega_bar_ms)+sigma_m1^2/2)/sigma_m1));

F_pis = normcdf((log(omega_bar_ms)+sigma_m1^2/2)/(sigma_m1));
F_pes = normcdf((log(omega_bar_es)+sigma_e1^2/2)/(sigma_e1));
%R_ms=R_tilde_Hs*(x_ms)/((Gamma_ms - mu_m*G_ms)*R_Hs);
%R_ms=R_mis;
R_tilde_Hs=R_ms*((Gamma_ms - mu_m*G_ms)*R_Hs)/(x_ms);
UL_m_1s=varphi_m*L_ms^(eta);
Lambda_ms=UL_m_1s/ws;
UC_m_1s=Lambda_ms;
C_ms=(1/UC_m_1s)/(1-hab);
%xi_ms=Lambda_ms/(rho_Hs*phi_H);
%xi_ms=(Lambda_ms-((1-F_pis)*betta_m*Lambda_ms*R_ms))/R_ms;
xi_ms=(Lambda_ms-((1-F_pis)*betta_m*Lambda_ms*R_ms))/((R_ms*(1-(1-F_pis)*(1-rp))));


%ZZHms =betta_m*v_m/(Lambda_ms*(1-betta_m*(1-Gamma_ms)*R_Hs)-xi_ms*((1-Gamma_Hs)*(Gamma_ms- mu_m*G_ms)*R_Hs));
%%ZZHms=betta_m*v_m/(Lambda_ms*(q_Hs) - betta_m*Lambda_ms*((1-G_ms)*R_Hs*q_Hs));
ZZHms=betta_m*v_m/(Lambda_ms*(q_Hs) - betta_m*Lambda_ms*(1-G_ms)*R_Hs*q_Hs-xi_ms*epsilonH*q_Hs*delta_H);

%-------------------------------------------------------------------------
%% SAVERS
%%-------------------------------------------------------------------------
L_ss=Ls-L_ms;

UL_s_1s=varphi_s*L_ss^(eta);

Lambda_ss=UL_s_1s/ws;

UC_s_1s = Lambda_ss;

C_ss=(1/UC_s_1s)/(1-hab);

Cs = C_ss + C_ms;

ZZHss =betta_s*v_s/(Lambda_ss*(1-(1-delta_H)*betta_s));

%q_Hs=(ZZHss+ZZHms)/Hs;
% Hs=1;
q_Hs = 1;
Hs=(ZZHss+ZZHms);

IHs = delta_H*Hs;

H_ss=ZZHss;
H_ms=ZZHms;

b_ms=x_ms*(H_ms*q_Hs)/R_ms;
%R_Ds=R_DDs;
%%---------------------------------------------------------------

%%-------------------------------------------------------------------------
%Tr_Hs = ((omega_bar_Hs - Gamma_Hs + mu_H*G_Hs)*R_tilde_Hs*(H_ms*q_Hs*x_ms)/R_ms);
% W_es = (1-Gamma_es)*R_Ks*q_Ks*Ks;%-a_e*Trs;
%         W_es = (1-Gamma_es)*R_Ks*q_Ks*Ks-a_e*Tr_Hs-a_e*((omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*R_tilde_Fs*(q_Ks*Ks-(1-chi_e)*W_es));
%VVe=(omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*R_tilde_Fs*(q_Ks*Ks);
%         W_es = (1-Gamma_es)*R_Ks*q_Ks*Ks-a_e*Tr_Hs-a_e*(VVe-(omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*R_tilde_Fs*(1-chi_e)*W_es);
%         W_es = (1-Gamma_es)*R_Ks*q_Ks*Ks-a_e*Tr_Hs-a_e*VVe+a_e*((omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*R_tilde_Fs*(1-chi_e)*W_es);
%         W_es*(1-a_e*((omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*R_tilde_Fs*(1-chi_e)) = (1-Gamma_es)*R_Ks*q_Ks*Ks-a_e*Tr_Hs-a_e*VVe;
%W_es=((1-Gamma_es)*R_Ks*q_Ks*Ks-a_e*Tr_Hs-a_e*VVe)/(1-a_e*((omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*R_tilde_Fs*(1-chi_e)));
%n_es = (1-chi_e)*W_es;
%b_es=q_Ks*Ks-(1-chi_e)*n_es;% b_es=q_Ks*Ks-(1-chi_e)*n_es;
%b_es=q_Ks*Ks-n_es;%
%Tr_Fs = ((omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*R_tilde_Fs*(q_Ks*Ks-(1-chi_e)*W_es));
W_es=((1-Gamma_es)*R_Ks*q_Ks*Ks);

n_es = (1-chi_e)*W_es;
b_es=q_Ks*Ks-n_es;
R_Fs= x_es*(q_Ks*Ks)/(b_es);
%xi_es=(R_Fs*q_Ks*(1-F_pes)-((((1-G_es)*R_Ks*q_Ks))*betta_s))/(epsilonF*q_Ks-R_Fs);
xi_es=(R_Fs*q_Ks*(1-F_pes)-((((1-G_es)*R_Ks*q_Ks))*betta_s))/(epsilonF*q_Ks*delta_K-((1-(1-F_pes)*(1-rpe))*R_Fs));
Tr_Hs = (omega_bar_Hs - Gamma_Hs + mu_H*G_Hs)*(((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms)));
Tr_Fs = (omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es)));
Trs = Tr_Fs + Tr_Hs;
%Trs = (omega_bar_Hs - Gamma_Hs + mu_H*G_Hs)*((((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms))+((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es))));
%W_bs = (phi_F*(q_Ks*Ks - (1-chi_e)*W_es) + phi_H*b_ms)/(1-chi_b);
%R_Ds=R_DDs;

Pr_Hs =(1-G_Hs)*(((1-F_pis)*b_ms*(R_ms)+G_ms*(1 - mu_m)*(b_ms*(R_ms)/omega_bar_ms)))+(1-G_Fs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*(R_Fs)/omega_bar_es)))-(1-F_pHs)*(b_ms*(1-phibs))*R_Ds-(1-F_pFs)*(b_es*(1-phibs))*R_Ds;

W_bs = Pr_Hs;
%rhos=(((1-G_Hs)*((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms)))+(1-G_Fs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es)))-(1-F_pHs)*(b_ms*(1-phi_H))*R_Ds-(1-F_pFs)*(b_es*(1-phi_F))*R_Ds)/(b_es*phi_F+b_ms*phi_H);
%W_bs = (phi_F*(q_Ks*Ks - (1-chi_e)*W_es) + phi_H*b_ms)*rhos;
%W_bs = (b_es*phi_F+b_ms*phi_H)*rhos;
C_bs=chi_b*W_bs;
%C_bs=W_bs-(b_es*phi_F+b_ms*phi_H);
%C_bs=chi_b*W_bs;
C_es=chi_e*W_es;
%C_bs=0;
%C_es=0;

%Trs =Tr_Fs + Tr_Hs;
%Trs =Tr_Fs + Tr_Hs+(C_bs+C_es);
%%-------------------------------------------------------------------------
E_Fs = phi_F*(q_Ks*Ks - n_es);


Ds = (ws*L_ss - delta_H*H_ss -C_ss + C_bs + C_es - a_s*Trs+PIs+PHs)/(1-R_DDs);
n_bs = (b_ms + (q_Ks*Ks - n_es)-Ds);%(rho_Fs*E_Fs + rho_Hs*(n_bs-E_Fs));

%W_bs =n_bs/(1-chi_b);

R_Fs= x_es*(q_Ks*Ks)/(q_Ks*Ks-n_es);
omega_bar_ms= x_ms/R_Hs; %(40) x_m
PHI_m = (1-Gamma_ms)*R_Hs*q_Hs*H_ms;


A1s=((betta_s*zeta1)*(Lambda_ss*((1-F_pHs)*(R_Ds)+nu*(phibs/phi)^(1-psib)/(b_ms+b_es))*((R_ms)^(tau))*b_ms))/(1-(betta_s*zeta1));


B1s=((betta_s*zeta1)*(Lambda_ss*((1-G_Hs)*(((1-F_pis))+(G_ms*(1 - mu_m)/omega_bar_ms)))*((((R_ms)^(tau))*b_ms))))/(1-(betta_s*zeta1));

C1s=((betta_s*zeta1)*(Lambda_ss*((1-F_pFs)*(R_Ds)+nu*(phibs/phi)^(1-psib)/(b_ms+b_es))*((R_Fs)^(tau))*b_es))/(1-(betta_s*zeta1));

D1s=(betta_s*zeta1)*(Lambda_ss*((1-G_Fs)*(((1-F_pes))+(G_es*(1 - mu_e)/omega_bar_es)))*((R_Fs)^(tau)*((b_es))))/(1-(betta_s*zeta1));

R_mis=(tau/(tau-1))*(A1s)/(B1s);
R_Fis=(tau/(tau-1))*(C1s)/(D1s);

R_Ds=R_DDs;
D_Fs = (1-phi_F)*E_Fs/phi_F;
D_Hs = (1-phi_H)*(n_bs-E_Fs)/phi_H;

av_defs = (D_Fs*def_rate_Fs + D_Hs*def_rate_Hs)/Ds;
%%-----------------------------------------------------------------
z(1) =R_mis-R_ms;%xi_bs/betta_s-rhos;%rhos-(((1-G_Hs)*((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms)))+(1-G_Fs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es)))-(1-F_pHs)*(b_ms*(1-phi_H))*R_Ds-(1-F_pFs)*(b_es*(1-phi_F))*R_Ds)/(b_es*phi_F+b_ms*phi_H);

z(2) = C_ms + q_Hs*H_ms - ws*L_ms - b_ms - PHI_m+a_b*Trs;
z(3) = -R_tilde_Fs +(Gamma_es-mu_e*G_es)*R_Ks*q_Ks*Ks/(q_Ks*Ks-n_es);%
z(4) = b_es-epsilonF*(Ks*delta_K)*q_Ks/((1-(1-F_pes)*(1-rpe))*R_Fs);%b_es*R_Fs-epsilonF*(Ks*q_Ks);%-xi_es + Gamma_e_1s/((1-Gamma_Fs)*(Gamma_e_1s- mu_e*G_e_1s));%(1-Gamma_es)*R_Ks + xi_es*((1-Gamma_Fs)*(Gamma_es- mu_e*G_es) *R_Ks - rho_Fs*phi_F);%(32)R_K
z(5) = b_ms-epsilonH*(H_ms*delta_H)*q_Hs/((1-(1-F_pis)*(1-rp))*R_ms);%b_ms*R_ms-epsilonH*(H_ms*q_Hs);%betta_m*Lambda_ms*Gamma_m_1s - xi_ms*(1-Gamma_Hs)*(Gamma_m_1s - mu_m*G_m_1s); %25
z(6)=n_bs-(1-chi_b)*W_bs;%n_bs-phi_F*(q_Ks*Ks - (1-chi_e)*W_es) + phi_H*b_ms;%rhos-(((1-G_Hs)*((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms)))+(1-G_Fs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es)))-(1-F_pHs)*(b_ms*(1-phi_H))*R_Ds-(1-F_pFs)*(b_es*(1-phi_F))*R_Ds)/(b_es*phi_F+b_ms*phi_H);%W_bs -n_bs/(1-chi_b);%rhos*(b_es*phi_F+b_ms*phi_H)-W_bs;%rho_Hs -(1-Gamma_Hs)*R_tilde_Hs/phi_H; %(40)R_tilde_Hs=rho_Hs*phi_H/(1-Gamma_Hs)
z(7)=omega_bar_Hs-((1-phi_H)*R_Ds)/R_tilde_Hs; %(38)
z(8)=R_Fis-R_Fs;%xi_bs-betta_s*((1-G_Fs)*(((1-F_pes)*(tau-1)*(R_Fs/infs1)+G_es*(1 - mu_e)*(tau-1)*((R_Fs/infs1)/omega_bar_es)))-(1-F_pFs)*((tau)*(1-phi_F))*R_Ds)/(phi_F*tau);
z(9)=omega_bar_Fs-((1-phi_F)*R_Ds)/R_tilde_Fs; %(37)
z(10)=phibs-n_bs/(b_ms+b_es);%R_mis-R_ms;%(betta_s*((1-G_Fs))*(((1-F_pes)+(1-mu_e)*G_es/omega_bar_es)*(1-taue))+betta_s*(((1-F_pFs)*((1-phi_F))*taue*R_Ds/R_Fs)+(phi_F*taue*rhos/R_Fs));%rhos-(((1-G_Hs)*((1-F_pis)*b_ms*(R_ms/infs1)+G_ms*(1 - mu_m)*(b_ms*(R_ms/infs1)/omega_bar_ms)))+(1-G_Fs)*(((1-F_pes)*b_es*R_Fs/infs1+G_es*(1 - mu_e)*(b_es*(R_Fs/infs1)/omega_bar_es)))-(1-F_pHs)*(b_ms*(1-phi_H))*R_Ds-(1-F_pFs)*(b_es*(1-phi_F))*R_Ds)/(b_es*phi_F+b_ms*phi_H);
z(11)=R_Ds-R_DDs/(1-pp*(av_defs));

result=solve(z,variables);
