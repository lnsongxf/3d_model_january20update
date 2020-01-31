
function z=Solve_SS_AdjRD(x,psib,nu,pu,zeta_m,zeta_F,def_rate_ss,tau_H,tau_F,rp,rpe,pp,hab,delta_H,a_e,a_s,a_b,q_Ks,As,PIs,PHs,R_Hs,R_DDs,alphaa,delta_K,betta_m,betta_s,phi_F,phi_H,phi,epsilonH,epsilonF,mu_m,mu_e,mu_B,sigma_e1,sigma_m1,sigma_B,varphi_s,varphi_m,v_s,v_m,chi_b,chi_e,eta,omikronH,omikronF);
error_flag=0;
R_Ks = x(1);
Ls   = x(2);
L_ms = x(3);
x_es = x(4);
x_ms = x(5);
omega_bar_Bs=x(6);
R_ms=x(7);
%omega_bar_Bs=xfs(8);
R_tilde_Fs=x(8);
phibs=x(9);

try
Gamma_Bs=normcdf((log(omega_bar_Bs) - sigma_B^2/2)/sigma_B) + omega_bar_Bs*(1-normcdf((log(omega_bar_Bs)+sigma_B^2/2)/sigma_B));
def_rate_Bs = def_rate_ss+normcdf((log(omega_bar_Bs)+sigma_B^2/2)/sigma_B)*400;
av_defs=def_rate_Bs;
R_Ds=R_DDs/(1-pp*(av_defs));

G_Bs   = normcdf((log(omega_bar_Bs)-sigma_B^2/2)/sigma_B);
F_pBs   = normcdf((log(omega_bar_Bs)+sigma_B^2/2)/sigma_B);

catch
   Gamma_Bs=Inf;
   def_rate_Bs=Inf;
   av_defs=def_rate_Bs;
   R_Ds=R_DDs/(1-pp*(av_defs));
   G_Bs=Inf;
   F_pBs=Inf;
   def_rate_Bs=Inf;
   error_flag=1;
end


%%-----------------------------------------------------------------
q_Hs=1;
r_Ks = R_Ks - (1-delta_K);
%rhos = 1/(1-chi_b);
YKs = r_Ks/alphaa;
Ks = ((YKs)/(Ls^(1-alphaa)))^(1/(alphaa-1));
Is = delta_K*Ks;

%%-------------------------------------------------------------------------
%% ENTREPRENEURS
%%-------------------------------------------------------------------------
Ys = As*Ks^(alphaa)*Ls^(1-alphaa);
ws = (1-alphaa)*Ys/Ls;
omega_bar_es = x_es/R_Ks;
try
G_es = normcdf((log(omega_bar_es)-sigma_e1^2/2)/sigma_e1 );
Gamma_es   = normcdf((log(omega_bar_es) - sigma_e1^2/2)/sigma_e1 ) + omega_bar_es*(1-normcdf((log(omega_bar_es)+sigma_e1^2/2)/sigma_e1 ));
catch
   
    G_es=Inf;
    Gamma_es=Inf;
    error_flag=1;
end

%xi_es=-(1-Gamma_es)*R_Ks/((1-Gamma_Bs)*(Gamma_es- mu_e*G_es) *R_Ks - rho_Fs*phi_F);
%%-------------------------------------------------------------------------
%% HOUSEHOLDS
%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------
%% BORROWERS
%%-------------------------------------------------------------------------
omega_bar_ms = x_ms/R_Hs;
try
omega_bar_ms = x_ms/R_Hs;
Gamma_ms   = normcdf((log(omega_bar_ms)- sigma_m1^2/2)/sigma_m1) + omega_bar_ms*(1-normcdf((log(omega_bar_ms)+sigma_m1^2/2)/sigma_m1));
G_ms   = normcdf((log(omega_bar_ms)-sigma_m1^2/2)/sigma_m1);
F_pis = normcdf((log(omega_bar_ms)+sigma_m1^2/2)/(sigma_m1));
F_pes = normcdf((log(omega_bar_es)+sigma_e1^2/2)/(sigma_e1));
catch
  Gamma_ms=Inf;
    G_ms=Inf;
F_pis=Inf;
    F_pes=Inf;
    error_flag=1;
end
%R_ms=R_tilde_Hs*(x_ms)/((Gamma_ms - mu_m*G_ms)*R_Hs);
%R_ms=R_mis;
R_tilde_Hs=R_ms*((Gamma_ms - mu_m*G_ms)*R_Hs)/(x_ms);
UL_m_1s=varphi_m*L_ms^(eta);
Lambda_ms=UL_m_1s/ws;
UC_m_1s=Lambda_ms;%10
C_ms=(1/UC_m_1s)/(1-hab);
%xi_ms=Lambda_ms/(rho_Hs*phi_H);
%xi_ms=(Lambda_ms-((1-F_pis)*betta_m*Lambda_ms*R_ms))/R_ms;

xi_ms=(Lambda_ms-((1-F_pis)*betta_m*Lambda_ms*(R_ms)))/(((R_ms)*(1-(1-F_pis)*(1-rp))));

ZZHms=betta_m*v_m/(Lambda_ms*(q_Hs) - betta_m*Lambda_ms*(1-G_ms)*R_Hs*q_Hs-xi_ms*epsilonH*q_Hs*delta_H);


%-------------------------------------------------------------------------
%% SAVERS
%%-------------------------------------------------------------------------
L_ss=Ls-L_ms;
UL_s_1s=varphi_s*L_ss^(eta);%11
Lambda_ss=UL_s_1s/ws;
UC_s_1s = Lambda_ss;%9UC
C_ss=(1/UC_s_1s)/(1-hab);
ZZHss =betta_s*v_s/(Lambda_ss*(1-(1-delta_H)*betta_s));

%q_Hs=(ZZHss+ZZHms)/Hs;
% Hs=1;
q_Hs = 1;
Hs=(ZZHss+ZZHms);
H_ss=ZZHss;
H_ms=ZZHms;
b_ms=x_ms*(H_ms*q_Hs)/(R_ms);
%R_Ds=R_DDs;
%%---------------------------------------------------------------

%%-------------------------------------------------------------------------

W_es=((1-Gamma_es)*R_Ks*q_Ks*Ks);

n_es = (1-chi_e)*W_es;
b_es=q_Ks*Ks-n_es;
R_Fs= (x_es*(q_Ks*Ks)/(b_es));

%Tr_Hs = (omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*(((1-F_pis)*b_ms*(R_ms)+G_ms*(1 - mu_m)*(b_ms*(R_ms)/omega_bar_ms)));
%Tr_Fs = (omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*(R_Fs)/omega_bar_es)));
Trs=(omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*((((1-F_pis)*b_ms*R_ms+(1-F_pes)*b_es*R_Fs+...
G_ms*(1 - mu_m)*(b_ms*(R_ms)/omega_bar_ms)+G_es*(1 - mu_e)*(b_es*(R_Fs)/omega_bar_es))));
%Trs = (omega_bar_Bs - Gamma_Bs + mu_B*G_Bs)*((((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms))+((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es))));
%W_bs = (phi_F*(q_Ks*Ks - (1-chi_e)*W_es) + phi_H*b_ms)/(1-chi_b);
%R_Ds=R_DDs;

% Pr_Hs =(1-G_Bs)*(((1-F_pis)*b_ms*(R_ms)+G_ms*(1 - mu_m)*(b_ms*(R_ms)/omega_bar_ms)))+(1-G_Bs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*(R_Fs)/omega_bar_es)))-(1-F_pBs)*(b_ms*(1-phibs))*R_Ds-(1-F_pBs)*(b_es*(1-phibs))*R_Ds;


%rhos=(((1-G_Bs)*((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms)))+(1-G_Bs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es)))-(1-F_pBs)*(b_ms*(1-phi_H))*R_Ds-(1-F_pBs)*(b_es*(1-phi_F))*R_Ds)/(b_es*phi_F+b_ms*phi_H);
%W_bs = (phi_F*(q_Ks*Ks - (1-chi_e)*W_es) + phi_H*b_ms)*rhos;
%W_bs = (b_es*phi_F+b_ms*phi_H)*rhos;

%C_bs=W_bs-(b_es*phi_F+b_ms*phi_H);
%C_bs=chi_b*W_bs;
C_es=chi_e*W_es;
%C_bs=0;
%C_es=0;

%Trs =Tr_Fs + Tr_Hs;
%Trs =Tr_Fs + Tr_Hs+(C_bs+C_es);
%%-------------------------------------------------------------------------




% n_bs = (b_ms + (q_Ks*Ks - n_es)-Ds);%(rho_Fs*E_Fs + rho_Hs*(n_bs-E_Fs));

%W_bs =n_bs/(1-chi_b);
R_Fs= (x_es*(q_Ks*Ks)/(q_Ks*Ks-n_es));
omega_bar_ms= x_ms/R_Hs;
PHI_m = (1-Gamma_ms)*R_Hs*q_Hs*H_ms;



% A1s=((betta_s*zeta1)*(Lambda_ss*((1-F_pBs)*(R_Ds)+nu*(phibs/phi)^(1-psib)/(b_ms+b_es))*((R_ms)^(tau))*b_ms))/(1-(betta_s*zeta1));

B1s=((betta_s*zeta_m)*(Lambda_ss*((1-G_Bs)*(((1-F_pis))+(G_ms*(1 - mu_m)/omega_bar_ms)))*((((R_ms)^(tau_H))*b_ms))))/(1-(betta_s*zeta_m));
%C1s=((betta_s*zeta1)*(Lambda_ss*((1-F_pBs)*(R_Ds)+nu*(phibs/phi)^(1-psib)/(b_ms+b_es))*((R_Fs)^(tau))*b_es))/(1-(betta_s*zeta1));

D1s=(betta_s*zeta_F)*(Lambda_ss*((1-G_Bs)*(((1-F_pes))+(G_es*(1 - mu_e)/omega_bar_es)))*((R_Fs)^(tau_F)*((b_es))))/(1-(betta_s*zeta_F));

% R_mis=(tau/(tau-1))*(A1s)/(B1s);
% R_Fis=(tau/(tau-1))*(C1s)/(D1s);

%R_Ds=R_DDs;


%  av_defs = (D_Fs*def_rate_Bs + D_Hs*def_rate_Bs)/Ds;

Pr_Hs = (1-G_Bs)*(((1-F_pis)*b_ms*(R_ms)+G_ms*(1 - mu_m)*(b_ms*(R_ms)/omega_bar_ms)))+(1-G_Bs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*(R_Fs)/omega_bar_es)))-(1-F_pBs)*(b_ms*(1-phibs))*R_Ds-(1-F_pBs)*(b_es*(1-phibs))*R_Ds;
W_bs= Pr_Hs;%14
C_bs	= chi_b*W_bs; %15
 Ds = (ws*L_ss - delta_H*H_ss -C_ss + C_bs + C_es - a_s*Trs+PIs+PHs)/(1-R_DDs); %16
 n_bs =(b_ms + (q_Ks*Ks - n_es)-Ds);
 A1s	 = ((betta_s*zeta_m)*(Lambda_ss*((1-F_pBs)*(R_Ds)+nu*(phibs/phi_H)^(1-psib)/(b_ms+b_es))*((R_ms)^(tau_H))*b_ms))/(1-(betta_s*zeta_m));
C1s      = ((betta_s*zeta_F)*(Lambda_ss*((1-F_pBs)*(R_Ds)+nu*(phibs/phi_F)^(1-psib)/(b_ms+b_es))*((R_Fs)^(tau_F))*b_es))/(1-(betta_s*zeta_F));
R_mis=(tau_H/(tau_H-1))*(A1s)/(B1s);
R_Fis=(tau_F/(tau_F-1))*(C1s)/(D1s);

%%-----------------------------------------------------------------
if error_flag==0
z(1) =R_mis-R_ms;%xi_bs/betta_s-rhos;%rhos-(((1-G_Bs)*((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms)))+(1-G_Bs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es)))-(1-F_pBs)*(b_ms*(1-phi_H))*R_Ds-(1-F_pBs)*(b_es*(1-phi_F))*R_Ds)/(b_es*phi_F+b_ms*phi_H);

z(2) = C_ms + q_Hs*H_ms - ws*L_ms - b_ms - PHI_m+a_b*Trs;
z(3) = -R_tilde_Fs +(Gamma_es-mu_e*G_es)*R_Ks*q_Ks*Ks/(q_Ks*Ks-n_es);%
z(4) = b_es-(epsilonF-omikronF*phi_F)*(Ks*delta_K)*q_Ks/((1-(1-F_pes)*(1-rpe))*R_Fs);%b_es*R_Fs-epsilonF*(Ks*q_Ks);%-xi_es + Gamma_e_1s/((1-Gamma_Bs)*(Gamma_e_1s- mu_e*G_e_1s));%(1-Gamma_es)*R_Ks + xi_es*((1-Gamma_Bs)*(Gamma_es- mu_e*G_es) *R_Ks - rho_Fs*phi_F);%(32)R_K
z(5) = b_ms-(epsilonH-omikronH*phi_H)*(H_ms*delta_H)*q_Hs/((1-(1-F_pis)*(1-rp))*R_ms);%b_ms*R_ms-epsilonH*(H_ms*q_Hs);%betta_m*Lambda_ms*Gamma_m_1s - xi_ms*(1-Gamma_Bs)*(Gamma_m_1s - mu_m*G_m_1s); %25
z(6)=n_bs-(1-chi_b)*W_bs;%n_bs-phi_F*(q_Ks*Ks - (1-chi_e)*W_es) + phi_H*b_ms;%rhos-(((1-G_Bs)*((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms)))+(1-G_Bs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es)))-(1-F_pBs)*(b_ms*(1-phi_H))*R_Ds-(1-F_pBs)*(b_es*(1-phi_F))*R_Ds)/(b_es*phi_F+b_ms*phi_H);%W_bs -n_bs/(1-chi_b);%rhos*(b_es*phi_F+b_ms*phi_H)-W_bs;%rho_Hs -(1-Gamma_Bs)*R_tilde_Hs/phi_H; %(40)R_tilde_Hs=rho_Hs*phi_H/(1-Gamma_Bs)
%z(7)=-omega_bar_Bs+(1-phibs)*R_Ds/((b_ms*R_tilde_Hs+b_es*R_tilde_Fs)/(b_es+b_ms)); %(38)
z(7)=-omega_bar_Bs+R_Ds*Ds/(R_tilde_Hs*b_ms+R_tilde_Fs*b_es);
z(8)=R_Fis-R_Fs;%xi_bs-betta_s*((1-G_Bs)*(((1-F_pes)*(tau-1)*(R_Fs/infs1)+G_es*(1 - mu_e)*(tau-1)*((R_Fs/infs1)/omega_bar_es)))-(1-F_pBs)*((tau)*(1-phi_F))*R_Ds)/(phi_F*tau);
z(9)=-phibs+n_bs/((((R_mis/R_ms)^(-tau_H))*b_ms+((R_Fis/R_Fs)^(-tau_F))*b_es));
 %z(9)=-phibs+n_bs/(b_ms+b_es);






else 
    z=Inf*ones(9,1);
end
% norm(z)

% if norm(imag(z))>0
%     z=Inf*ones(10,1);
% end

% z=sum(z.^2);
end

%;

