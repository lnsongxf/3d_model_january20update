//%%%%%**********************************************************************
//% Codes to solve the 3D MODEL
//% UPON USE PLEASE CITE: "Capital Regulation in a Macroeconomic Model with Three Layers of Default"
//% by  Laurent Clerc, Alexis Derviz, Caterina Mendicino, Stephane Moyen, Kalin Nikolov, 
//% Livio Stracca, Javier Suarez, Alexandros P. Vardoulakis, 
//% International Journal of Central Banking, June 2015, Pages 9-63.
//%%%%%**********************************************************************
//% THIS PROGRAM WAS TESTED WITH (AND REQUIRES KNOWLEDGE OF)
//% 1) DYNARE 4.4.3 (http://www.cepremap.cnrs.fr/dynare/)
//% 2) MATLAB R2012b (http://www.mathworks.com/)
//%%%%%**********************************************************************
//% Contact:
//% CATERINA MENDICINO: caterina.mendicino1@ecb.int
//% KALIN NIKOLOV     : kalin.Nikolov@ecb.int
//% DOMINIK SUPERA    : dominik.supera@gmail.com
//%%%%%**********************************************************************
//% Latest update: 22/02/2016 by DOMINIK SUPERA
//%%%%%**********************************************************************


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//1. Preamble 
//The preamble consists of the some declarations to setup the endogenous
//and exogenous variables, the parameters and assign values to these parameters.
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//***************************************************
// DECLARATION ENDOGENOUS VARIABLES IN THE MODEL
//***************************************************

var 
b_e        // entrepreneurial debt
g_H        // housing investment adjustment cost
g_H_1      // first derivative of g_H with respect to I_H
IH         // housing investment 
PH         // Profit of housing capital producing firm
A          // productivity
b_m        // mortgage debt
C          // aggregate consumption
C_m        // consumption of borrowers (impatient households)
C_s        // consumption of savers (patient households)
D          // aggregate deposits
def_rate_e // corporate default rate
def_rate_m // mortgage default rate//def_rate_F // default rate of corporate banks
def_rate_H // default rate of mortgage banks//E_F        // equity invested in corporate banks
def_rate_F
epsilonH
epsilonF
G_e        // share of entrepreneurial capital belonging to firms that default (BGG parameter)
G_e_1      // First derivative of G_e with respect to omega_e (BGG parameter)//G_F        // share of corporate loans belonging to corporate banks that default (BGG parameter)
G_F        // share of mortgage loans belonging to mortgage banks that default (BGG parameter)
F_pF        //cdf for bank default
G_H        // share of mortgage loans belonging to mortgage banks that default (BGG parameter)
F_pH
g_I        // Capital investment adjustment cost 
g_I_1      // First derivative of g_I with respect to I_K
G_m        // share of housing belonging to households that default (BGG parameter)
G_m_1      // First derivative of G_m with respect to omega_m (BGG parameter)
Gamma_e    // Share of gross corporate revenues going to the bank (BGG parameter)
Gamma_e_1  // First derivative of Gamma_e with respect to omega_e(BGG parameter)//Gamma_F    // Share of gross corporate bank revenues going to depositors(BGG parameter)//Gamma_F_1  // First derivative of Gamma_F with respect to omega_F(BGG parameter)
Gamma_F
Gamma_F_1
Gamma_H    // Share of gross mortgage bank revenues going to depositors(BGG parameter)
Gamma_H_1  // First derivative of Gamma_H with respect to omega_H(BGG parameter)
Gamma_m    // Share of gross returns from housing going to the bank(BGG parameter)
Gamma_m_1  // First derivative of Gamma_m with respect to omega_m(BGG parameter)
F_pi       // cdf of probability of default
F_pe        //probability of defaultng households
//F_pi_1
//F_pe_1
H          // Aggregate housing supply
H_m        // Housing used by borrowers
H_s        // Housing used by savers
I          // Aggregate corporate investment
K          // Aggregate capital stock 
L          // Aggregate labour supply
L_m        // Borrowers' labour supply
L_s        // Savers' labour supply       
Lambda_m   // LM on the budget constraint in the household borrower's problem
Lambda_s   // LM on the budget constraint in the household saver's problem
n_b        // Bankers' net worth
n_e        // Entrepreneurs' net worth
omega_bar_e	// Idiosyncratic productivity shock below which the corporate borrower defaults
omega_bar_F // Idiosyncratic corporate bank loan return shock below which the corporate bank defaults	
omega_bar_H	// Idiosyncratic mortgage bank loan return shock below which the mortgage bank defaults
omega_bar_m // Idiosyncratic housing return shock below which the mortgage borrower defaults	
PI          // Profit of business capital producing firm 
q_H         // Housing price
q_K         // Capital price
R_D         // Deposit interest rate
R_F         // Corporate loan interest rate NOMINAL
R_H         // Aggregate financial return on housing (i.e. excluding imputed rents)
r_K         // Capital rental rate
R_K         // Capital rate of return
R_m         // Mortgage interest rate NOMINAL
inf1         
R_tilde_F 	// Aggregate return on a diversified corporate loan portfolio (i.e. portfolio return after accounting for loan losses)
R_tilde_H 	// Aggregate return on a diversified housing loan portfolio (i.e. portfolio return after accounting for loan losses)
//rho_F       // Corporate bank return on equity
//rho_H       // Mortgage bank return on equity
Tr_H
Tr_F
Tr          // Total deposit insurance payments to banks//Tr_F        // deposit insurance payments to corporate banksTr_H        // deposit insurance payments to mortgage banks
UC_m        // Borrowers' utility from consumption
UC_m_1      // Borrowers' marginal utility from consumption
UC_s        // Savers' utility from consumption
UC_s_1      // Savers' marginal utility from consumption
UH_m        // Borrowers' utility from housing services
UH_m_1      // Borrowers' marginal utility from housing services
UH_s        // Savers' utility from housing services
UH_s_1      // Savers' marginal utility from housing services
UL_m        // Borrowers' disutility from work
UL_m_1      // Borrowers' marginal disutility from work
UL_s        // Savers' disutility from work
UL_s_1      // Savers' marginal disutility from work
Util_m      // Total borrower utility
Util_s      // Total saver utility
w           // Wage rate
W_b         // Bankers' wealth (pre-dividend)
Pr_H        //Profit of the bank
//rho         //return on bak equity
W_e         // Entrepreneurs' wealth (pre-dividend)
x_e         // Corporate leverage
x_m         // Household leverage
xi_e        // LM on bank's participation constraint in the entrepreneurs' problem
xi_m        // LM on bank's participation constraint in the household borrower's problem
//xi_b          // LM on bank's BS constraint 
Y           // Output
Y_obs       // All _obs variables are transformations to aid plotting. Transformations may differ across variables depending on what makes sense.
R_D_obs     // Please look at particular _obs variable definition in the code to see the exact transformation.
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
//E_F_obs
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
//rho_F_obs
//rho_H_obs
Tr_obs
//Tr_F_obs
//Tr_H_obs
w_obs
x_e_obs
x_m_obs
Vs           // Value function of household savers
Vm           // Value function of household borrowers
Ve           // Value function of entrepreneurs
Vb           // Value function of bankers
EJ           // Shock (housing preference)
EK           // Shock (Capital Investment)
EH           // Shock (Housing Investment)
ESe          // Shock (entrepreneur risk)
ESm          // Shock (housing risk)
ESH          // Shock (mortgage bank risk)
ESF          // Shock (corporate bank risk)
EWe          // Shock (Entrepreneur net-worth)
EWb          // Shock (Banker net-worth)
EdH          // Shock (housing depreciation)
EdK          // Shock (capital depreciation)
//ERW        //RW shock
phi_F        // Minimum capital ratio for corporate banks
phi_H        // Minimum capital ratio for mortgage banks
phi
b_e_obs
bsp_F        // Corporate loan return spread
bsp_H        // Household loan return spread
Y_net        // Output net of default costs
Y_net_obs
Y_net_2
Y_net_2_obs
res_chk
deltaH       // Housing depreciation rate
deltaK       // Capital depreciation rate
GDP          // GDP
GDP_obs
H_obs
IH_obs
res_H
m_e          // Verification cost of entrepreneurs
m_m          // Verification cost of borrowing households
av_def       // Average default for all banks
RDsp         // Deposit rate spread over the risk free rate
R_DD         // Deposit rate return net of bank default costs
W_m          // Household borrower net worth
W_m_obs
b_tot        // Total debt in the economy
b_tot_obs
C_b          // Dividend payments by bankers
C_e          // Dividend payments by entrepreneurs
A1
B1
C1
D1
R_mi
R_Fi
phib
;

parameters sigma_epsiHd sigma_epsiHk sigma_epsiA sigma_epsiJ sigma_epsiK sigma_epsiH 
sigma_epsiSe sigma_epsiSm sigma_epsiSF sigma_epsiSH sigma_epsiWb sigma_epsiRW sigma_epsiWe pp  
hab Cyphi_H Cyphi_F phi_Fs phi_Hs alphaa delta_K delta_H betta_m betta_s mu_m mu_e mu_F mu_H 
sigma_e1 sigma_m1 sigma_F sigma_H varphi_s varphi_m v_s v_m chi_b chi_e eta  a_e a_s a_b   psi_i psi_h 
rhoA rhoJ rhoK rhoSe rhoSm rhoSF rhoSH rhoWb rhoWe rhoHd rhoHk rhoRW zeta1 zetae epsilonH1s epsilonF1s tau taue rp rpe phiinf1 kappa nu psib phis;
    
varexo epsiA  epsiJ epsiK epsiSe epsiSm epsiSF epsiWb  epsiWe epsiH epsiHd ELTV ;
//epsiRW
//************************************************************
// Parametrization
//************************************************************

load par_ThreeDAdj_RD_mod  ;
filename = 'par_ThreeDAdj_RD_mod.mat';
m = matfile(filename,'Writable',true);

m.pp=0;
m.phis=0.2;
m.psib=10;
m.nu=0.3;
m.chi_b=0.1;
m.chi_e=0.1;
m.alphaa=0.2;
m.rp=0.035;
m.rpe=0.04;
m.delta_H=0.035;
m.delta_K=0.035;
m.tau=50;
m.taue=50;
m.zeta1=0.45;
m.zetae=0.45;
m.kappa=0;
m.betta_m=0.97;
m.sigma_epsiRW=25;
m.epsilonH1s=0.8;
m.epsilonF1s=0.8;
m.rhoRW=0.97;
m.phi_Fs = 0.2;
m.phi_Hs = 0.2;
m.v_m=0.25;
m.v_s=0.25;
m.mu_H=0.3;
m.mu_m=0.25;
m.mu_e=0.25;
m.sigma_m1=0.12;
m.sigma_e1=0.2;
m.sigma_epsiRW=25;
m.phiinf1=1.01;


set_param_value('alphaa',alphaa);
set_param_value('delta_K',delta_K);
set_param_value('delta_H',delta_H);
set_param_value('betta_m',betta_m);
set_param_value('betta_s',betta_s);
set_param_value('mu_m',mu_m);
set_param_value('mu_e',mu_e);
set_param_value('mu_F',mu_F);
set_param_value('mu_H',mu_H);
set_param_value('sigma_e1',sigma_e1);
set_param_value('sigma_m1',sigma_m1);
set_param_value('sigma_F',sigma_F);
set_param_value('sigma_H',sigma_H);
set_param_value('varphi_s',varphi_s);
set_param_value('varphi_m',varphi_m);
set_param_value('v_s',v_s);
set_param_value('v_m',v_m);
set_param_value('chi_b',chi_b);
set_param_value('chi_e',chi_e);
set_param_value('eta',eta);
set_param_value('a_e',a_e);
set_param_value('a_s',a_s);
set_param_value('a_b',a_b);
set_param_value('psi_i',psi_i);
set_param_value('psi_h',psi_h);
set_param_value('rhoA',rhoA);
set_param_value('rhoJ',rhoJ);
set_param_value('rhoA',rhoA);
set_param_value('rhoK',rhoK);
set_param_value('rhoSe',rhoSe);
set_param_value('rhoSm',rhoSm);
set_param_value('rhoSF',rhoSF);
set_param_value('rhoSH',rhoSH);
set_param_value('rhoWe',rhoWe);
set_param_value('rhoWb',rhoWb);
set_param_value('rhoHd',rhoHd);
set_param_value('hab',hab);
set_param_value('rhoHk',rhoHk);
set_param_value('phi_Fs',phi_Fs);
set_param_value('phi_Hs',phi_Hs);
set_param_value('Cyphi_H',Cyphi_H);
set_param_value('Cyphi_F',Cyphi_F);
set_param_value('pp',pp);
set_param_value('sigma_epsiA',sigma_epsiA);
set_param_value('sigma_epsiJ',sigma_epsiJ);
set_param_value('sigma_epsiK',sigma_epsiK);
set_param_value('sigma_epsiH',sigma_epsiH);
set_param_value('sigma_epsiSe',sigma_epsiSe);
set_param_value('sigma_epsiSm',sigma_epsiSm);
set_param_value('sigma_epsiSF',sigma_epsiSF);
set_param_value('sigma_epsiSH',sigma_epsiSH);
set_param_value('sigma_epsiWb',sigma_epsiWb);
set_param_value('sigma_epsiWe',sigma_epsiWe);
set_param_value('sigma_epsiHd',sigma_epsiHd);
set_param_value('sigma_epsiHk',sigma_epsiHk);
set_param_value('sigma_epsiRW',sigma_epsiRW);
set_param_value('rhoRW',rhoRW);
set_param_value('zeta1',zeta1);
set_param_value('zetae',zetae);
set_param_value('tau',tau);
set_param_value('taue',taue);
set_param_value('rp',rp);
set_param_value('rpe',rpe);
set_param_value('epsilonH1s',epsilonH1s);
set_param_value('epsilonF1s',epsilonF1s);
set_param_value('kappa',kappa);
set_param_value('phiinf1',phiinf1);
set_param_value('psib',psib);
set_param_value('phis',phis);
set_param_value('nu',nu);
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// 2. DECLARATION OF THE MODEL
// It starts with the instruction "model;" and ends with "end;", in between all equilibrium
// conditions are written exactly the way we write it “by hand”. 
// Note that:
// - if x is decided in period t then we simply write x. 
// - when the variable is decided in t-1, we write x(-1).
// - when a variable is decided in the next period, t + 1, we write x(+1).
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//************************************************************
model;
//************************************************************
// Households 
//************************************************************
//*******************************
// Utility functions
//*******************************
//1
UC_s = log(C_s-hab*C_s(-1)); 
//2
UC_m = log(C_m-hab*C_m(-1)); 
//3
UL_s = varphi_s*L_s^(1+eta)/(1+eta); 
//4
UL_m = varphi_m*L_m^(1+eta)/(1+eta); 
//5
UH_s = EJ*v_s*log(H_s(-1)); 
//6
UH_m = EJ*v_m*log(H_m(-1)); 
//7
Util_s = UC_s - UL_s + UH_s;
//8
Util_m = UC_m - UL_m + UH_m;
//9
UC_s_1 = 1/(C_s-hab*C_s(-1));
//10
UC_m_1 = 1/(C_m-hab*C_m(-1));
//11
UL_s_1 = varphi_s*L_s^(eta);
//12
UL_m_1 = varphi_m*L_m^(eta);
//13
UH_s_1 = EJ*v_s/(H_s(-1));
//14
UH_m_1 = EJ*v_m/(H_m(-1));
//*******************************
// WELFARE
//*******************************
// 15
Vs=Util_s+betta_s*Vs(+1);
// 16
Vm=Util_m+betta_m*Vm(+1);
// 17
Ve=n_e;
// 18
Vb=n_b;
//*******************************
// Savers
//*******************************
// 19 foc C_s
Lambda_s = UC_s_1;
// 20 foc L_s
UL_s_1 = w*Lambda_s;
// 21 foc d
Lambda_s = betta_s*Lambda_s(1)*R_DD(1);
//22 Definition of R_DD (deposit return net of bank default costs)
R_DD=R_D(-1)/inf1*(1-pp*(av_def));
// 23 foc H_s
Lambda_s*(q_H) = betta_s*UH_s_1(1) + betta_s*Lambda_s(1)*(1-deltaH(1))*q_H(1);
// 24 Budget Constraint
  C_s - C_e - C_b + q_H*(H_s-(1-deltaH)*H_s(-1))+ D = w*L_s +R_DD*D(-1) - Tr*a_s + PI + PH +(1-EWe)*(1-Gamma_e)*R_K*q_K(-1)*K(-1) + (1-EWb)*(W_b); 
//*******************************
// Borrowers
//*******************************
//25 Default cut off borrowers depends on leverage
omega_bar_m = x_m(-1)/R_H;
//26 Household leverage definition
x_m = (R_m/inf1(+1))*b_m/(H_m*q_H);
//27 Housing rate of return
R_H = (1-deltaH)*q_H/q_H(-1);
//28 foc C_m    
Lambda_m = UC_m_1;
//29 foc L_m
UL_m_1 = w*Lambda_m;
// 30 foc X_m
//- betta_m*Lambda_m(1)*Gamma_m_1(1) + xi_m*((1-Gamma_H(1))*(Gamma_m_1(1) - m_m*G_m_1(1)))=0;
//- betta_m*Lambda_m(1)*Gamma_m_1(1) + xi_m*((1-Gamma_H(1))*(Gamma_m_1(1) - m_m*G_m_1(1)))-xi_m*zeta1*(R_m-R_m(-1))+xi_m(+1)*zeta1*(R_m(+1)-R_m)=0;
//- betta_m*Lambda_m(1)*Gamma_m_1(1)*b_m + xi_m*((1-Gamma_H(1))*(Gamma_m_1(1) - m_m*G_m_1(1)))*b_m -xi_m*zeta1*(R_m-steady_state(R_m))=0;
//- betta_m*Lambda_m(1)*Gamma_m_1(1) + xi_m*((1-Gamma_H(1))*(Gamma_m_1(1) - m_m*G_m_1(1)))=0;

// 30
//b_m*R_m=epsilonH*H_m*q_H;
(b_m-b_m(-1)*(1-rp)*(1-F_pi))*R_m/inf1(+1)=epsilonH*(H_m-H_m(-1)*(1-delta_H))*q_H(+1);
// 31
//b_e*R_F=epsilonF*K*q_K;
(b_e-b_e(-1)*(1-rpe)*(1-F_pe))*R_F/inf1(+1)=epsilonF*(K-K(-1)*(1-delta_K))*q_K(+1);
// 32
epsilonH=epsilonH1s*ELTV;
// 33
epsilonF=epsilonF1s*ELTV;
// 34 foc H_m
//betta_m*UH_m_1(1) - Lambda_m*(q_H) + betta_m*Lambda_m(1)*(1-Gamma_m(1))*R_H(1)*q_H +xi_m*((1-Gamma_H(1))*(Gamma_m(1)- m_m*G_m(1))*R_H(1))*q_H= 0; 
//betta_m*UH_m_1(1) - Lambda_m*(q_H) + betta_m*Lambda_m(1)*(1-Gamma_m(1))*R_H(1)*q_H +xi_m*((1-Gamma_H(1))*(Gamma_m(1)- m_m*G_m(1))*R_H(1))*q_H-xi_m*zeta1*(R_m-R_m(-1))*R_m/H_m+xi_m(+1)*zeta1*(R_m(+1)-R_m)*R_m/H_m= 0;                                         
//betta_m*UH_m_1(1) - Lambda_m*(q_H) + betta_m*Lambda_m(1)*(1-Gamma_m(1))*R_H(1)*q_H +xi_m*((1-Gamma_H(1))*(Gamma_m(1)- m_m*G_m(1))*R_H(1))*q_H-xi_m*zeta1*((R_m-steady_state(R_m)))*R_m/H_m+xi_b*(epsilon*(q_H))= 0;                                         
//betta_m*UH_m_1(1) - Lambda_m*(q_H) + betta_m*Lambda_m(1)*(1-Gamma_m(1))*R_H(1)*q_H +xi_m*((1-Gamma_H(1))*(Gamma_m(1)- m_m*G_m(1))*R_H(1))*q_H-xi_m*zeta1*((R_m-steady_state(R_m)))*R_m/H_m= 0;                                         
//betta_m*UH_m_1(1) - Lambda_m*(q_H) + betta_m*Lambda_m(1)*(1-Gamma_m(1))*R_H(1)*q_H +xi_m*((1-Gamma_H(1))*(Gamma_m(1)- m_m*G_m(1))*R_H(1))*q_H= 0;  //needs to be changed                                       
//betta_m*UH_m_1(1) - Lambda_m*(q_H) + betta_m*Lambda_m(1)*(1-G_m(1))*R_H(1)*q_H = 0;  //needs to be changed                                       

betta_m*UH_m_1(1) - Lambda_m*(q_H) + betta_m*Lambda_m(1)*((1-G_m(1))*R_H(1)*q_H)+(xi_m)*(epsilonH*(q_H(+1)))-xi_m(+1)*(epsilonH*(q_H(+2)))*(1-deltaH) = 0;  //needs to be changed                                       
//betta_m*UH_m_1(1) - Lambda_m*(q_H) + betta_m*Lambda_m(1)*(1-Gamma_m(1))*R_H(1)*q_H +xi_m*((1-Gamma_H(1))*(Gamma_m(1)- m_m*G_m(1))*R_H(1))*q_H+(betta_m*Lambda_m(1)*Gamma_m_1(1) - xi_m*((1-Gamma_H(1))*(Gamma_m_1(1) - m_m*G_m_1(1))))*(b_m*R_m)/H_m+xi_m*(epsilonH*(q_H))-xi_m(+1)*(epsilonH*(q_H(+1)))*(1-deltaH)= 0;                                         

//betta_m*v_m/ZZHm-(Lambda_ms*(q_Hs) + betta_m*Lambda_ms*((1-G_ms)*R_Hs*q_Hs+xi_ms*epsilonH*q_Hs))s;



// 35 foc b_m 
//Lambda_m-xi_m*rho_H(1)*phi_H=0;
//Lambda_m-xi_m*rho_H(1)*phi_H+xi_m*zeta1*(R_m-R_m(-1))*R_m/b_m-xi_m(+1)*zeta1*(R_m(+1)-R_m)*R_m/b_m(R_m-steady_state(R_m))=0;
//Lambda_m-xi_m*rho_H(1)*phi_H+xi_m*zeta1*((R_m-steady_state(R_m)))*R_m/b_m-xi_b*R_m=0;
//Lambda_m-xi_m*rho_H(1)*phi_H+xi_m*zeta1*((R_m-steady_state(R_m)))*R_m/b_m=0;
//Lambda_m-xi_m*rho_H(1)*phi_H+xi_m*zeta1*((R_m-steady_state(R_m)))*R_m=0;///needs to be changed
Lambda_m-betta_m*Lambda_m(1)*(R_m/inf1(+1)*(1-F_pi(+1)))-xi_m*R_m/inf1(+1)+xi_m(+1)*(1-rp)*(1-F_pi(+1))*R_m(+1)/inf1(+2)=0;///needs to be changed


// 36 Budget Constraint Borrowers
C_m + q_H*H_m - (1-Gamma_m)*R_H*q_H(-1)*H_m(-1) = w*L_m + b_m-Tr*a_b;
//************************************************************
// Entrepreneurs
//************************************************************
// 37 Rate of return to capital
R_K = (r_K+(1-deltaK)*q_K)/q_K(-1);
// 38 Default cut off entrepreneurs depends on leverage
omega_bar_e = x_e(-1)/R_K;
// 39 Firm leverage definition
x_e = (R_F/inf1(+1))*(q_K*K-n_e)/(q_K*K);
//  FOC x_e
//- Gamma_e_1(1) + xi_e*((1-Gamma_F(1))*(Gamma_e_1(1)- m_e*G_e_1(1)));             
// 40 Foc k
//(1-Gamma_e(1))*R_K(1) + xi_e*((1-Gamma_F(1))*(Gamma_e(1)- m_e*G_e(1)) *R_K(1) - rho_F(1)*phi_F);//need to be changed for pc
(((1-G_e(1))*R_K(1)*q_K))*betta_s*(Lambda_s(+1)/Lambda_s)-(q_K*R_F/inf1(+1)*(1-F_pe(+1)))+xi_e*(epsilonF*q_K(+1)-R_F/inf1(+1))-xi_e(+1)*(epsilonF*q_K(+2)*(1-deltaK)-(1-rpe)*(1-F_pe(+1))*R_F(+1)/inf1(+2))=0 ;//need to be changed for pc
//xi_es=(R_Fs*q_Ks*(1-F_pes)-((((1-G_es)*R_Ks*q_Ks))*betta_s))/(epsilonF*q_Ks*delta_K-(R_Fs*(1-(1-F_pes)*(1-rpe))));

//(1-Gamma_e(1))*R_K(1)=0;//need to be changed for pc


// 41 Wealth dynamics entrepreneurs  
W_e= EWe*(1-Gamma_e)*R_K*q_K(-1)*K(-1);
// 42 Evolution individual net worth 
n_e = (1-chi_e)*W_e;
// 43 Corporate dividends
C_e=chi_e*W_e;


//***********************************************************
// Bankers
//************************************************************
// 42 Equalization of expected rates of return on equity invested in corporate and mortgage banks
//rho_F(1) = rho_H(1);
// 44 Wealth of bankers before dividends
//W_b = (rho_F*E_F(-1) + rho_H*(n_b(-1)-E_F(-1)));
//W_b(1)=rho*(b_e*phi_F+b_m*phi_H);
//W_b=rho*(b_e(-1)*phi_F(-1)+b_m(-1)*phi_H(-1));
W_b=Pr_H;

//Pr_H =(1-Gamma_H)*R_tilde_H*b_m(-1)+(1-Gamma_H)*R_tilde_F*(q_K(-1)*K(-1)-n_e(-1));
//Pr_H =(1-G_H)*(((1-F_pi)*b_m(-1)*R_m(-1)+G_m*(1 - mu_m)*(b_m(-1)*R_m(-1)/omega_bar_m)))+(1-G_F)*(((1-F_pe)*b_e(-1)*R_F(-1)+G_e*(1 - mu_e)*(b_e(-1)*R_F(-1)/omega_bar_e)))-(1-F_pH)*(b_m(-1)*(1-phi_H))*R_D(-1)-(1-F_pF)*(b_e(-1)*(1-phi_F))*R_D(-1);


// 45 Banker net worth allocated to lending
//n_b = (1-chi_b)*n_b(-1)+(1-ghi_b)*(W_b-n_b(-1));
n_b = (1-chi_b)*W_b;
// 45 Market clearing in corporate bank equity market
//E_F = phi_F*(q_K*K - n_e); 
//n_b=phi_H*b_m+phi_F*b_e;

//xi_b=rho(+1)*betta_s;
// 46 Bank dividends
C_b=chi_b*(W_b);
//************************************************************
// Banks
//************************************************************
// Definitions
// 47 Default threshold corporate bank
omega_bar_F = (1-phi_F(-1))*(R_D(-1)/inf1)/R_tilde_F;
// 48 Default threshold mortgage bank
omega_bar_H = (1-phi_H(-1))*(R_D(-1)/inf1)/R_tilde_H;
//omega_bar_H = (D(-1)*R_D(-1))/(R_tilde_H*b_m(-1)+R_tilde_H*b_e(-1));


// 49 Rate of return on corporate bank equity
//rho_F = (1-Gamma_F)*R_tilde_F/phi_F(-1);
// 50 Rate of return on mortgage bank equity
//rho_H =(1-Gamma_H)*R_tilde_H/phi_H(-1);

// 49
//Pr_H(+1) =(1-Gamma_H)*R_tilde_H*b_m(-1)+(1-Gamma_H)*R_tilde_F*(q_K(-1)*K(-1)-n_e(-1));//Profit equation
//Pr_H(+1) =((1-G_H(+1))*((1-F_pi(+1))*b_m*R_m/inf1(+1)+G_m(+1)*(1 - mu_m)*(b_m*(R_m/inf1(+1))/omega_bar_m(+1))))+(1-G_F(+1))*(((1-F_pe(+1))*b_e*R_F/inf1(+1)+G_e(+1)*(1 - mu_e)*(b_e*(R_F/inf1(+1))/omega_bar_e)))-(1-F_pH(+1))*(b_m*(1-phi_H))*R_D/inf1(+1)-(1-F_pF(+1))*(b_e*(1-phi_F))*R_D/inf1(+1)-rho(+1)*(b_e*phi_F+b_m*phi_H);
Pr_H(+1) =((1-G_H(+1))*((1-F_pi(+1))*b_m*R_m/inf1(+1)+G_m(+1)*(1 - mu_m)*(b_m*(R_m/inf1(+1))/omega_bar_m(+1))))+(1-G_F(+1))*(((1-F_pe(+1))*b_e*R_F/inf1(+1)+G_e(+1)*(1 - mu_e)*(b_e*(R_F/inf1(+1))/omega_bar_e)))-(1-F_pH(+1))*(b_m*(1-phib))*R_D/inf1(+1)-(1-F_pF(+1))*(b_e*(1-phib))*R_D/inf1(+1);

//Pr_Hs =(1-G_Hs)*(((1-F_pis)*b_ms*(R_ms/infs1)+G_ms*(1 - mu_m)*(b_ms*(R_ms/infs1)/omega_bar_ms)))+(1-G_Fs)*(((1-F_pes)*b_es*R_Fs/infs1+G_es*(1 - mu_e)*(b_es*(R_Fs/infs1)/omega_bar_es)))-(1-F_pHs)*(b_ms*(1-phi_H))*R_Ds-(1-F_pFs)*(b_es*(1-phi_F))*R_Ds-rhos*(b_es*phi_F+b_ms*phi_H)-0.054733173756168;


//R_tilde_H = (Gamma_m - mu_m*G_m)*R_H*q_H(-1)*H_m(-1)/b_m(-1);
//rho_H =(1-Gamma_H)*R_tilde_H/phi_H(-1);
//Pr_H =((1-G_H)*((1-F_pi)*b_m(-1)*R_m(-1)+G_m*(1 - mu_m)*(b_m(-1)*R_m(-1)/omega_bar_m)))+(1-G_F)*(((1-F_pe)*b_e(-1)*R_F(-1)+G_e*(1 - mu_e)*(b_e(-1)*R_F(-1)/omega_bar_e)))-(1-F_pH)*(b_m(-1)*(1-phi_H(-1)))*R_D(-1)-(1-F_pF)*(b_e(-1)*(1-phi_F(-1)))*R_D(-1)-rho(-1)*(b_e(-1)*phi_F(-1)+b_m(-1)*phi_H(-1));
//rho=(((1-G_H)*((1-F_pi)*b_m(-1)*R_m(-1)+G_m*(1 - mu_m)*(b_m(-1)*R_m(-1)/omega_bar_m)))+(1-G_F)*(((1-F_pe)*b_e(-1)*R_F(-1)+G_e*(1 - mu_e)*(b_e(-1)*R_F(-1)/omega_bar_e)))-(1-F_pH)*(b_m(-1)*(1-phi_H(-1)))*R_D(-1)-(1-F_pF)*(b_e(-1)*(1-phi_F(-1)))*R_D(-1))/(n_b(-1));

//Pr_H(+1)=0;

//Pr_Hs(+1) =(1-G_Hs)*(((1-F_pis)*b_ms*R_ms+G_ms*(1 - mu_m)*(b_ms*R_ms/omega_bar_ms(+1))))+(1-G_Fs)*(((1-F_pes)*b_es*R_Fs+G_es*(1 - mu_e)*(b_es*R_Fs/omega_bar_es)))-(1-F_pHs)*(b_ms*(1-phi_H))*R_Ds-(1-F_pFs)*(b_es*(1-phi_F))*R_Ds;
//rho=(((1-G_H(+1))*((1-F_pi(+1))*b_m*R_m+G_m(+1)*(1 - mu_m)*(b_m*R_m/omega_bar_m(+1))))+(1-G_F(+1))*(((1-F_pe(+1))*b_e*R_F+G_e(+1)*(1 - mu_e)*(b_e*R_F/omega_bar_e)))-(1-F_pH(+1))*(b_m*(1-phi_H))*R_D-(1-F_pF(+1))*(b_e*(1-phi_F))*R_D)/(b_e*phi_F+b_m*phi_H);
//rho=(((1-G_H)*((1-F_pi)*b_m(-1)*R_m(-1)+G_m*(1 - mu_m)*(b_m(-1)*R_m(-1)/omega_bar_m)))+(1-G_F)*(((1-F_pe)*b_e(-1)*R_F(-1)+G_e*(1 - mu_e)*(b_e(-1)*R_F(-1)/omega_bar_e)))-(1-F_pH)*(b_m(-1)*(1-phi_H(-1)))*R_D(-1)-(1-F_pF)*(b_e(-1)*(1-phi_F(-1)))*R_D(-1))/(n_b(-1));

///50 FOC for bank mortgage lending
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_H(+1))*((1-F_pi(+1)*R_m+G_m(+1)*(1 - mu_m)*(R_m/omega_bar_m(+1)))))-(1-F_pH(+1))*((1-phi_H))*R_D)-xi_b*phi_H=0;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_H(+1))*(((1-F_pi(+1))*(1-tau)+G_m(+1)*(1 - mu_m)*(((1-tau))/omega_bar_m(+1)))))+(1-F_pH(+1))*((1-phi_H))*(R_D/R_m)*tau+(phi_H*tau*rho/R_m))+xi_b*phi_H*tau/R_m-zeta1*((R_m-R_m(-1))/R_m(-1))*(R_m/R_m(-1))+betta_s*(Lambda_s(+1)/Lambda_s)*zeta1*((R_m(+1)-R_m))*((R_m(+1)/R_m)^2)*(b_m(+1)/b_m)=0;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_H(+1))*(((1-F_pi(+1))*(1-tau)+G_m(+1)*(1 - mu_m)*(((1-tau))/omega_bar_m(+1)))))+(1-F_pH(+1))*((1-phi_H))*(R_D/R_m)*tau+(phi_H*tau*rho/R_m))+xi_b*phi_H*tau/R_m-zeta1*(R_m-steady_state(R_m))=0;//+betta_s*(Lambda_s(+1)/Lambda_s)*zeta1*((R_m(+1)-R_m))*((R_m(+1)/R_m)^2)*(b_m(+1)/b_m)=0;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_H(+1))*(((1-F_pi(+1))*(1-tau)+G_m(+1)*(1 - mu_m)*(((1-tau))/omega_bar_m(+1)))))+(1-F_pH(+1))*((1-phi_H))*(R_D/R_m)*tau+(phi_H*tau*rho/R_m))+xi_b*phi_H*tau/R_m-zeta1*(R_m-steady_state(R_m))=0;//+betta_s*(Lambda_s(+1)/Lambda_s)*zeta1*((R_m(+1)-R_m))*((R_m(+1)/R_m)^2)*(b_m(+1)/b_m)=0;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_H(+1))*((1-F_pi(+1)*R_m+G_m(+1)*(1 - mu_m)*(R_m/omega_bar_m(+1)))))-(1-F_pH(+1))*((1-phi_H))*R_D)-xi_b*phi_H=0;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_H(+1))*(((1-F_pi(+1))*(1-tau)+G_m(+1)*(1 - mu_m)*(((1-tau))/omega_bar_m(+1)))))+(1-F_pH(+1))*((1-phi_H))*(R_D/R_m)*tau+(phi_H*tau*rho(+1)/R_m))+xi_b*phi_H*tau/R_m-zeta1*((R_m-R_m(-1)))*(R_m)+betta_s*(Lambda_s(+1)/Lambda_s)*zeta1*((R_m(+1)-R_m))*(b_m(+1)/b_m)*(R_m(+1))=0;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_H(+1))*(((1-F_pi(+1))*(1-tau)+G_m(+1)*(1 - mu_m)*(((1-tau))/omega_bar_m(+1)))))+(1-F_pH(+1))*((1-phi_H))*((R_D/inf1(+1))/(R_m/inf1(+1)))*tau+(phi_H*tau*rho(+1)/(R_m/inf1(+1))))+xi_b*phi_H*tau/(R_m/inf1(+1))-zeta1*((R_m-R_m(-1)))+betta_s*(Lambda_s(+1)/Lambda_s)*zeta1*((R_m(+1)-R_m))=0;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_H(+1))*(((1-F_pi(+1))*(1-tau)+G_m(+1)*(1 - mu_m)*(((1-tau))/omega_bar_m(+1)))))+(1-F_pH(+1))*((1-phi_H))*((R_D/inf1(+1))/(R_m/inf1(+1)))*tau+(phi_H*tau*rho(+1)/(R_m/inf1(+1))))-zeta1*((R_m-R_m(-1)))+betta_s*(Lambda_s(+1)/Lambda_s)*zeta1*((R_m(+1)-R_m))=0;
//(1-G_Fs)*(((1-F_pes)*(tau-1)*(R_Fs/infs1)+G_es*(1 - mu_e)*(tau-1)*((R_Fs/infs1)/omega_bar_es)))-(1-F_pFs)*((tau)*(1-phi_F))*R_Ds-((1-G_Hs)*(((1-F_pis)*(tau-1)*(R_ms/infs1)+G_ms*(1 - mu_m)*(tau-1)*((R_ms/infs1)/omega_bar_ms)))-(1-F_pHs)*((tau)*(1-phi_H))*R_Ds);

//R_m=1.03;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_H(+1))*(((1-F_pi(+1))*(1-tau)+G_m(+1)*(1 - mu_m)*(((1-tau))/omega_bar_m(+1)))))+(1-F_pH(+1))*((1-phi_H))*(R_D/R_m)*tau+(phi_H*tau*rho(+1)/R_m))+xi_b*phi_H*tau/R_m-zeta1*((R_m-steady_state(R_m)));//+betta_s*(Lambda_s(+1)/Lambda_s)*zeta1*((R_m(+1)-R_m))*(b_m(+1)/b_m)*(R_m(+1))=0;

//Pr_H(+1) =((1-G_H(+1))*((1-F_pi(+1))*b_m*R_m/inf1(+1)+G_m(+1)*(1 - mu_m)*(b_m*(R_m/inf1(+1))/omega_bar_m(+1))))+(1-G_F(+1))*(((1-F_pe(+1))*b_e*R_F/inf1(+1)+G_e(+1)*(1 - mu_e)*(b_e*(R_F/inf1(+1))/omega_bar_e)))-(1-F_pH(+1))*(b_m*(1-phi_H))*R_D/inf1(+1)-(1-F_pF(+1))*(b_e*(1-phi_F))*R_D/inf1(+1)-rho(+1)*(b_e*phi_F+b_m*phi_H);
//Pr_H(+1) =((1-G_H(+1))*((1-F_pi(+1))*(((R_mi/R_m)^(-tau))*b_m)*R_m/inf1(+1)+G_m(+1)*(1 - mu_m)*((((R_mi/R_m)^(-tau))*b_m)*(R_m/inf1(+1))/omega_bar_m(+1))))+(1-G_F(+1))*(((1-F_pe(+1))*(((R_Fi/R_F)^(-tau))*b_e)*R_F/inf1(+1)+G_e(+1)*(1 - mu_e)*((((R_Fi/R_F)^(-tau))*b_e)*(R_F/inf1(+1))/omega_bar_e)))-(1-F_pH(+1))*((((R_mi/R_m)^(-tau))*b_m)*(1-phi_H))*R_D/inf1(+1)-(1-F_pF(+1))*((((R_Fi/R_F)^(-tau))*b_e)*(1-phi_F))*R_D/inf1(+1)-rho(+1)*((((R_Fi/R_F)^(-tau))*b_e)*phi_F+(((R_mi/R_m)^(-tau))*b_m)*phi_H);

//R_mi=(tau/(tau-1))*((betta_s*zeta1)^s*(Lambda_s(+1)*((1-F_pH(+1))*((1-phi_H))*(R_D/inf1(+1))+rho(+1)*phi_H)*((R_m)^(tau))*b_m))/((betta_s*zeta1)^s*Lambda_s(+1)*((1-G_H(+1))*(((1-F_pi(+1)))+(G_m(+1)*(1 - mu_m)/omega_bar_m(+1))))*((((R_m)^(tau))*b_m)/inf1(+1)));
//R_Fi=(tau/(tau-1))*((betta_s*zeta1)^s*(Lambda_s(+1)*((1-F_pH(+1))*((1-phi_F))*(R_D/inf1(+1))+rho(+1)*phi_F)*((R_F)^(tau))*b_e))/((betta_s*zeta1)^s*Lambda_s(+1)*((1-G_F(+1))*(((1-F_pe(+1)))+(G_e(+1)*(1 - mu_m)/omega_bar_e(+1))))*((((R_F)^(tau))*b_e)/inf1(+1)));

phib=n_b/((((R_mi/R_m)^(-tau))*b_m+((R_Fi/R_F)^(-tau))*b_e));

//A=(betta_s*zeta1)^s*(Lambda_s(+1)*((1-F_pH(+1))*((1-phi_H))*(R_D/inf1(+1))+rho(+1)*phi_H)*((R_m)^(tau))*b_m);
A1 = (betta_s*zeta1)*(Lambda_s(+1)*((1-F_pH(+1))*(R_D/inf1(+1))+nu*(phib/phi)^(1-psib)/(((R_mi/R_m)^(-tau))*b_m+((R_Fi/R_F)^(-tau))*b_e))*((R_m)^(tau))*b_m)+(betta_s*zeta1)*A1(+1);
//nu*(phib/phi)^(1-psib)/(((R_mi/R_m)^(-tau))*b_m+((R_Fi/R_F)^(-tau))*b_e)
//A1 = (betta_s)*(Lambda_s(+1)*((1-F_pH(+1))*((1-phi_H))*(R_D/inf1(+1))+xi_b*phi_H/betta_s)*((R_m)^(tau))*b_m)/inf1(+1);
//A1s=((betta_s*zeta1)*(Lambda_ss*((1-F_pHs)*((1-phi_H))*(R_Ds/infs1)+xi_bs*phi_H/betta_s)*((R_ms)^(tau))*b_ms))/(1-(betta_s*zeta1));
//A1s=(((1-F_pHs)*((1-phi_H))*(R_Ds/infs1)+xi_bs*phi_H/betta_s)*((R_ms)^(tau))*b_ms);
//B1s=(((1-G_Hs)*(((1-F_pis))+(G_ms*(1 - mu_m)/omega_bar_ms)))*((((R_ms)^(tau))*b_ms)/infs1));

//xi_bs=betta_s*(((1-G_Hs)*(((1-F_pis)*(tau-1)*(R_ms/infs1)+G_ms*(1 - mu_m)*(tau-1)*((R_ms/infs1)/omega_bar_ms)))-(1-F_pHs)*((tau)*(1-phi_H))*R_Ds))/(phi_H*tau);



//B=(betta_s*zeta1)^s*Lambda_s(+1)*((1-G_H(+1))*(((1-F_pi(+1)))+(G_m(+1)*(1 - mu_m)/omega_bar_m(+1))))*((((R_m)^(tau))*b_m)/inf1(+1));
B1=(betta_s*zeta1)*Lambda_s(+1)*((1-G_H(+1))*(((1-F_pi(+1)))+(G_m(+1)*(1 - mu_m)/omega_bar_m(+1))))*((((R_m)^(tau))*b_m)/inf1(+1))+(betta_s*zeta1)*B1(+1);
//B1=(betta_s)*Lambda_s(+1)*((1-G_H(+1))*(((1-F_pi(+1)))+(G_m(+1)*(1 - mu_m)/omega_bar_m(+1))))*((((R_m)^(tau))*b_m)/inf1(+1));
//B1s=((betta_s*zeta1)*(Lambda_ss*((1-G_Hs)*(((1-F_pis))+(G_ms*(1 - mu_m)/omega_bar_ms)))*((((R_ms)^(tau))*b_ms)/infs1)))/(1-(betta_s*zeta1));

C1=(betta_s*zeta1)*(Lambda_s(+1)*((1-F_pF(+1))*(R_D/inf1(+1))+(nu*(phib/phi)^(1-psib)/(((R_mi/R_m)^(-tau))*b_m+((R_Fi/R_F)^(-tau))*b_e)))*((R_F)^(tau))*b_e)+(betta_s*zeta1)*C1(+1);
//C1=(betta_s)*(Lambda_s(+1)*((1-F_pH(+1))*((1-phi_F))*(R_D/inf1(+1))+xi_b*phi_F/betta_s)*((R_F)^(tau))*b_e)/inf1(+1);


D1=(betta_s*zeta1)*Lambda_s(+1)*((1-G_F(+1))*(((1-F_pe(+1)))+(G_e(+1)*(1 - mu_e)/omega_bar_e(+1))))*((((R_F)^(tau))*b_e)/inf1(+1))+(betta_s*zeta1)*D1(+1);
//D1=(betta_s)*Lambda_s(+1)*((1-G_F(+1))*(((1-F_pe(+1)))+(G_e(+1)*(1 - mu_e)/omega_bar_e(+1))))*((((R_F)^(tau))*b_e)/inf1(+1));

R_mi=(tau/(tau-1))*(A1)/(B1);
R_Fi=(tau/(tau-1))*(C1)/(D1);
R_m=((1-zeta1)*R_mi^(1-tau)+zeta1*R_m(-1)^(1-tau))^(1/(1-tau));
R_F=((1-zeta1)*R_Fi^(1-tau)+zeta1*R_F(-1)^(1-tau))^(1/(1-tau));


//b_mi=((R_mi/R_m)^(-tau))*b_m
//b_ei=((R_Fi/R_F)^(-tau))*b_e

//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_H(+1))*(((1-F_pi(+1))*(1)+G_m(+1)*(1 - mu_m)*(((1))/omega_bar_m(+1)))))-(1-F_pH(+1))*((1-phi_H))*(R_D/R_m)*1-(phi_H*1*rho(+1)/R_m))-xi_b*phi_H*1/R_m;//-zeta1*((R_m-R_m(-1)))*(R_m)+betta_s*(Lambda_s(+1)/Lambda_s)*zeta1*((R_m(+1)-R_m))*(b_m(+1)/b_m)*(R_m(+1))=0;
//xi_bs=betta_s*(((((1-G_Hs)*(((1-F_pis)*(1-tau*0)+G_ms*(1 - mu_m)*((1-tau*0)/omega_bar_ms))))))-((R_Ds/R_ms)*(1)*((1-F_pHs)*(1-phi_H))-(phi_H*1*rhos/R_ms)))*R_ms/(phi_H*1);


//51FOC for bank business lending
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_F(+1))*(((1-F_pe(+1))*R_F+G_e(+1)*(1 - mu_e)*(R_F/omega_bar_e(+1)))))-(1-F_pF(+1))*((1-phi_F))*R_D)-xi_b*phi_F=0;

//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_F(+1))*(((1-F_pe(+1))*(1-taue)+G_e(+1)*(1 - mu_e)*((1-taue)/omega_bar_e(+1)))))+(1-F_pF(+1))*((1-phi_F))*(R_D/R_F)*taue+(phi_F*taue*rho/R_F))+xi_b*phi_F*taue/R_F-zetae*((R_F-R_F(-1))/R_F(-1))*(R_F/R_F(-1))+betta_s*(Lambda_s(+1)/Lambda_s)*zetae*((R_F(+1)-R_F))*((R_F(+1)/R_F)^2)*(b_e(+1)/b_e)=0;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_F(+1))*(((1-F_pe(+1))*(1-taue)+G_e(+1)*(1 - mu_e)*((1-taue)/omega_bar_e(+1)))))+(1-F_pF(+1))*((1-phi_F))*(R_D/R_F)*taue+(phi_F*taue*rho/R_F))+xi_b*phi_F*taue/R_F-zetae*(R_F-steady_state(R_F))=0;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_F(+1))*(((1-F_pe(+1))*(1-taue)+G_e(+1)*(1 - mu_e)*((1-taue)/omega_bar_e(+1)))))+(1-F_pF(+1))*((1-phi_F))*(R_D/R_F)*taue+(phi_F*taue*rho(+1)/R_F))+xi_b*phi_F*taue/R_F-zetae*((R_F-R_F(-1)))+betta_s*(Lambda_s(+1)/Lambda_s)*zetae*((R_F(+1)-R_F))*(b_e(+1)/b_e)*R_F(+1)=0;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_F(+1))*(((1-F_pe(+1))*(1-taue)+G_e(+1)*(1 - mu_e)*((1-taue)/omega_bar_e(+1)))))+(1-F_pF(+1))*((1-phi_F))*(R_D/R_F)*taue+(phi_F*taue*rho(+1)/R_F))+xi_b*phi_F*taue/R_F-zetae*((R_F-steady_state(R_F)));//+betta_s*(Lambda_s(+1)/Lambda_s)*zetae*((R_F(+1)-R_F))*(b_e(+1)/b_e)*R_F(+1)=0;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_F(+1))*(((1-F_pe(+1))*(1-taue)+G_e(+1)*(1 - mu_e)*((1-taue)/omega_bar_e(+1)))))+(1-F_pF(+1))*((1-phi_F))*((R_D/inf1(+1))/(R_F/inf1(+1)))*taue+(phi_F*taue*rho(+1)/(R_F/inf1(+1))))+xi_b*phi_F*taue/(R_F/inf1(+1))-zetae*((R_F(+1)-R_F(-1)))+betta_s*(Lambda_s(+1)/Lambda_s)*zetae*((R_F(+1)-R_F))=0;
//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_F(+1))*(((1-F_pe(+1))*(1-taue)+G_e(+1)*(1 - mu_e)*((1-taue)/omega_bar_e(+1)))))+(1-F_pF(+1))*((1-phi_F))*((R_D/inf1(+1))/(R_F/inf1(+1)))*taue+(phi_F*taue*rho(+1)/(R_F/inf1(+1))))-zetae*((R_F(+1)-R_F(-1)))+betta_s*(Lambda_s(+1)/Lambda_s)*zetae*((R_F(+1)-R_F))=0;


//betta_s*(Lambda_s(+1)/Lambda_s)*(((1-G_F(+1))*(((1-F_pe(+1))*(1)+G_e(+1)*(1 - mu_e)*((1)/omega_bar_e(+1)))))-(1-F_pF(+1))*((1-phi_F))*(R_D/R_F)*1-(phi_F*1*rho(+1)/R_F))-xi_b*phi_F*1/R_F;//-zetae*((R_F-R_F(-1)))+betta_s*(Lambda_s(+1)/Lambda_s)*zetae*((R_F(+1)-R_F))*(b_e(+1)/b_e)*R_F(+1)=0;
//betta_s*((1-G_Fs))*(((1-F_pes)+(1-mu_e)*G_es/omega_bar_es)*(1))-betta_s*((1-F_pFs)*((1-phi_F))*1*R_Ds/R_Fs)-(phi_F*1*rhos/R_Fs)-xi_bs*phi_F*1/R_Fs
//52
//n_b=b_m*phi_H+b_e*phi_F;
//betta_s*(Lambda_s(+1)/Lambda_s)*(1-Gamma_H(+1))*R_tilde_F(+1)=xi_bk;
//xi_bk=betta_s*(Lambda_s(+1)/Lambda_s)*R_DD(1);

// 53 Rate of return on household mortgage loans
R_tilde_H = (Gamma_m - mu_m*G_m)*R_H*q_H(-1)*H_m(-1)/b_m(-1);
// 54 Rate of return on corporate loans
R_tilde_F = (Gamma_e-mu_e*G_e)*R_K*q_K(-1)*K(-1)/(q_K(-1)*K(-1)-n_e(-1));
// 55 Balance Sheet Bank
n_b + D = b_m + (q_K*K - n_e); 
//************************************************************
// Consumption good production
//************************************************************
// 56 Output
Y = A*K(-1)^(alphaa)*L^(1-alphaa);
// 57 Capital rental rate
r_K = alphaa*Y/K(-1);
// 58 Wage rate
w = (1-alphaa)*Y/L;
//************************************************************
// Capital good production with CEE adjustment costs 
// (depend on I/I(-1)) - Christiano, Eichenbaum and Evans (1995) - JPE
//************************************************************
// 59 Capital adj cost function
g_I = EK*psi_i/2*(I/I(-1)-1)^2;
// 60 foc of (57) wrt I
g_I_1 = EK*psi_i*(I/I(-1)-1);
// 61 Capital stock evolution
K = (1-deltaK)*K(-1) + I*(1-g_I );
// 62 Capital price
q_K = 1 + g_I + I/I(-1)*g_I_1 - betta_s*Lambda_s(1)/Lambda_s*(I(1)/I)^2*g_I_1(1);
// 63 Flow profit of capital producers
PI = q_K*I - (1 + g_I)*I; 
//*******************************
// Inside equity market
//*******************************
// 64 Market clearing in bank equity market
//(1-chi_b)*W_b = phi_F*(q_K*K - (1-chi_e)*W_e) + phi_H*b_m; 
//*******************************
// Goods market
//*******************************
// 65 Aggregate consumption definition
C = C_s + C_m;
//*******************************
// Labour market 
//*******************************
// 66 Aggregate labour supply definition
L = L_s + L_m;
//*******************************
// Housing supply with CEE Adj costs
//*******************************
// 67 Capital adj cost function
g_H = EH*psi_h/2*(IH/IH(-1)-1)^2;
// 68 foc of (63) wrt IH
g_H_1 = EH*psi_h*(IH/IH(-1)-1);
// 69 Housing stock evolution
H = (1-deltaH)*H(-1) + IH*(1-g_H);
// 70 Housing price with CEE Adj costs
q_H = 1 + g_H + IH/IH(-1)*g_H_1 - betta_s*Lambda_s(1)/Lambda_s*(IH(1)/IH)^2*g_H_1(1);
// 71 Flow profits by housing producers 
PH= q_H*IH - (1 + g_H)*IH; 
// 72 Aggregate housing market clearing (Supply (LHS) = Demand (RHS))
H = H_m + H_s; 
//*******************************
// Deposit insurance agency
//*******************************
//  Deposit insurance transfers to corporate bank depositors whose banks have gone bankrupt
//Tr_F = (omega_bar_F - Gamma_F + mu_F*G_F)*R_tilde_F*(q_K(-1)*K(-1)-(1-chi_e)*W_e(-1));
// 70 Deposit insurance transfers to household bank depositors whose banks have gone bankrupt
Tr_H = (omega_bar_H - Gamma_H + mu_H*G_H)*((((1-F_pi)*b_m(-1)*R_m(-1)/inf1+G_m*(1 - mu_m)*(b_m(-1)*(R_m(-1)/inf1)/omega_bar_m))));
Tr_F = (omega_bar_F - Gamma_F + mu_F*G_F)*(((1-F_pe)*b_e(-1)*R_F(-1)/inf1+G_e*(1 - mu_e)*(b_e(-1)*(R_F(-1)/inf1)/omega_bar_e)));

//Tr = (omega_bar_H - Gamma_H + mu_H*G_H)*(R_tilde_H*b_m+R_tilde_F*b_e);
//Tr = (omega_bar_H - Gamma_H + mu_H*G_H)*(R_tilde_H*b_m+R_tilde_F*b_e);
// 73
//Tr = (omega_bar_H - Gamma_H + mu_H*G_H)*((((1-F_pi)*b_m(-1)*R_m(-1)+G_m*(1 - mu_m)*(b_m(-1)*R_m(-1)/omega_bar_m))+((1-F_pe)*b_e(-1)*R_F(-1)+G_e*(1 - mu_e)*(b_e(-1)*R_F(-1)/omega_bar_e))));
// 71 Total deposit insurance transfers
Tr = Tr_F + Tr_H;



//**************************************************************
//Taylor rule
//**************************************************************

(R_D/steady_state(R_D))=((R_D(-1)/steady_state(R_D))^kappa)*((inf1/steady_state(inf1))^phiinf1)^(1-kappa);




//************************************************************
// Distribution functions
//************************************************************
// 74
Gamma_m   = normcdf((log(omega_bar_m)- (ESm*sigma_m1)^2/2)/(ESm*sigma_m1)) + omega_bar_m*(1-normcdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1)));
// 75
Gamma_e   = normcdf((log(omega_bar_e) - (ESe*sigma_e1)^2/2)/(ESe*sigma_e1)) + omega_bar_e*(1-normcdf((log(omega_bar_e)+ESe*sigma_e1 ^2/2)/(ESe*sigma_e1)));
// 76
Gamma_F   = normcdf((log(omega_bar_F) - (ESF*sigma_F)^2/2)/(ESF*sigma_F)) + omega_bar_F*(1-normcdf((log(omega_bar_F)+ESF*sigma_F^2/2)/(ESF*sigma_F)));
// 77
Gamma_H   = normcdf((log(omega_bar_H)-(ESH*sigma_H)^2/2)/(ESH*sigma_H)) + omega_bar_H*(1-normcdf((log(omega_bar_H)+ESH*sigma_H^2/2)/(ESH*sigma_H)));
// 78
Gamma_m_1 = (normpdf((log(omega_bar_m)-(ESm*sigma_m1)^2/2)/(ESm*sigma_m1))-omega_bar_m*normpdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1)))/(ESm*sigma_m1*omega_bar_m) + (1-normcdf((log(omega_bar_m)+(ESm*sigma_m1)^2/2)/(ESm*sigma_m1)));
// 79
Gamma_e_1 = (normpdf((log(omega_bar_e)-(ESe*sigma_e1) ^2/2)/(ESe*sigma_e1))-omega_bar_e*normpdf((log(omega_bar_e)+ESe*sigma_e1 ^2/2)/(ESe*sigma_e1)))/(ESe*sigma_e1 *omega_bar_e) + (1-normcdf((log(omega_bar_e)+(ESe*sigma_e1)^2/2)/(ESe*sigma_e1)));
// 80
Gamma_F_1 = (normpdf((log(omega_bar_F)-(ESF*sigma_F)^2/2)/(ESF*sigma_F))-omega_bar_F*normpdf((log(omega_bar_F)+ESF*sigma_F^2/2)/(ESF*sigma_F)))/(ESF*sigma_F*omega_bar_F) + (1-normcdf((log(omega_bar_F)+(ESF*sigma_F)^2/2)/(ESF*sigma_F)));
// 81
Gamma_H_1 = (normpdf((log(omega_bar_H)-(ESH*sigma_H)^2/2)/(ESH*sigma_H))-omega_bar_H*normpdf((log(omega_bar_H)+ESH*sigma_H^2/2)/(ESH*sigma_H)))/(ESH*sigma_H*omega_bar_H) + (1-normcdf((log(omega_bar_H)+(ESH*sigma_H)^2/2)/(ESH*sigma_H)));
// 82
G_m   = normcdf((log(omega_bar_m)-ESm*sigma_m1^2/2)/(ESm*sigma_m1));
// 83
G_e   = normcdf((log(omega_bar_e)-ESe*sigma_e1 ^2/2)/(ESe*sigma_e1));
// 84
G_F   = normcdf((log(omega_bar_F)-ESF*sigma_F^2/2)/(ESF*sigma_F));
// 85
G_H   = normcdf((log(omega_bar_H)-ESH*sigma_H^2/2)/(ESH*sigma_H));
// 86
G_m_1 = normpdf((log(omega_bar_m)-ESm*sigma_m1^2/2)/(ESm*sigma_m1))/(ESm*sigma_m1*omega_bar_m);
// 87
F_pi = normcdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1));//probability of defaultng households
//F_pi_1 = normpdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1))/(ESm*sigma_m1*omega_bar_m);
// 88
F_pe = normcdf((log(omega_bar_e)+ESm*sigma_e1^2/2)/(ESm*sigma_e1));//probability of defaultng households
//F_pe_1 = normpdf((log(omega_bar_e)+ESm*sigma_e1^2/2)/(ESm*sigma_e1))/(ESe*sigma_e1 *omega_bar_e);//probability of defaultng households
// 89
F_pH   = normcdf((log(omega_bar_H)+sigma_H^2/2)/sigma_H);
// 90
F_pF   = normcdf((log(omega_bar_F)+sigma_F^2/2)/sigma_F);
// 91
G_e_1 = normpdf((log(omega_bar_e)-ESe*sigma_e1 ^2/2)/(ESe*sigma_e1))/(ESe*sigma_e1 *omega_bar_e);

// 92 Capital requirement on mortgage banks
//exp(phi_H-phi_Hs) = (b_m/steady_state(b_m))^(Cyphi_H);
phi_H = phi_Hs;

// 93 Capital requirement on corporate banks
//exp(phi_F-phi_Fs) = (b_e/steady_state(b_e))^(Cyphi_F);
phi_F = phi_Fs;
phi = phis;

// Capital depreciation
// 94
deltaK = delta_K+EdK;
// 95 Housing depreciation
deltaH = delta_H+EdH;

//*********************
//****Other definitions
//*********************
// NFC loans and total loans
//96
b_e=(q_K*K-n_e);
//97
b_tot=b_e+b_m;
// Return on loan portfolio spreads
//98
bsp_H = 400*(R_tilde_H(1)-R_D);
//99
bsp_F = 400*(R_tilde_F(1)-R_D);
//GDP DEFINITIONS (Y_net is the main definition of GDP used in the model)
//100
Y_net = C + (chi_b)*W_b + (chi_e)*W_e + I + IH;
//101
Y_net_2 = Y - (R_D(-1)*pp*av_def*D(-1)+m_e*G_e*R_K*q_K(-1)*K(-1) + m_m*G_m*R_H*q_H(-1)*H_m(-1)+ mu_F*G_H*(Gamma_e - m_e*G_e)*R_K*q_K(-1)*K(-1) + mu_H*G_H*(Gamma_m - m_m*G_m)*R_H*q_H(-1)*H_m(-1));
GDP = Y + UL_s_1*H_s/Lambda_s + UL_m_1*H_m/Lambda_m;
// Monitoring cost
m_e = mu_e;
m_m = mu_m;
//Average default of banks
av_def = ((1-phi_F)*b_e*def_rate_F + (1-phi_H)*b_m*def_rate_H)/D;

// Deposit spread
RDsp = 400*(R_D-R_DD(+1));
// Borrowers wealth
W_m = (1-Gamma_m)*R_H*q_H(-1)*H_m(-1) + w*L_m;

//************************************************************
// Reporting adjunct variables
//************************************************************
def_rate_m = normcdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1))*400;
def_rate_e = normcdf((log(omega_bar_e)+ESe*sigma_e1^2/2)/(ESe*sigma_e1))*400;
def_rate_H = normcdf((log(omega_bar_H)+ESH*sigma_H^2/2)/(ESH*sigma_H))*400;
def_rate_F = normcdf((log(omega_bar_F)+ESF*sigma_F^2/2)/(ESF*sigma_F))*400;
Y_obs = log(Y/steady_state(Y))*100; 
R_D_obs = (R_D-steady_state(R_D))*400; 
R_m_obs = (R_m-steady_state(R_m))*400;
R_H_obs = log(R_H/steady_state(R_H))*400;
R_F_obs = (R_F-steady_state(R_F))*400;
H_m_obs = log(H_m/steady_state(H_m))*100;
H_s_obs = log(H_s/steady_state(H_s))*100;
b_m_obs = log(b_m/steady_state(b_m))*100;
b_e_obs = log(b_e/steady_state(b_e))*100;
C_obs = log(C/steady_state(C))*100;
C_m_obs = log(C_m/steady_state(C_m))*100;
C_s_obs = log(C_s/steady_state(C_s))*100;
D_obs = log(D/steady_state(D))*100;
//E_F_obs = log(E_F/steady_state(E_F))*100;
I_obs = log(I/steady_state(I))*100;
K_obs = log(K/steady_state(K))*100;
L_obs = log(L/steady_state(L))*100;
L_m_obs = log(L_m/steady_state(L_m))*100;
L_s_obs = log(L_s/steady_state(L_s))*100;
n_b_obs = log(n_b/steady_state(n_b))*100;
n_e_obs = log(n_e/steady_state(n_e))*100;
q_H_obs = log(q_H/steady_state(q_H))*100;
q_K_obs = log(q_K/steady_state(q_K))*100;
r_K_obs = log(r_K/steady_state(r_K))*400;
R_K_obs = (R_K-steady_state(R_K))*400;
R_tilde_F_obs = (R_tilde_F-steady_state(R_tilde_F))*400;
R_tilde_H_obs = (R_tilde_H-steady_state(R_tilde_H))*400;
//rho_F_obs = (rho_F-steady_state(rho_F))*400;
//rho_H_obs = (rho_H-steady_state(rho_H))*400;
Tr_obs = log(Tr/steady_state(Tr))*100;
//Tr_F_obs = log(Tr_F/steady_state(Tr_F))*100;
//Tr_H_obs = log(Tr_H/steady_state(Tr_H))*100;
w_obs = log(w/steady_state(w))*100;
x_e_obs = log(x_e/steady_state(x_e))*100;
x_m_obs = log(x_m/steady_state(x_m))*100;
b_tot_obs=log(b_tot/steady_state(b_tot))*100;
W_m_obs = log(W_m/steady_state(W_m))*100;
GDP_obs = log(GDP/steady_state(GDP))*100;
H_obs = log(H/steady_state(H))*100;
IH_obs = log(IH/steady_state(IH))*100;
Y_net_2_obs = log(Y_net_2/steady_state(Y_net_2))*100;
Y_net_obs = log(Y_net/steady_state(Y_net))*100;


// Residuals
res_chk = (R_D(-1)*pp*av_def*D(-1)+C + chi_b*W_b + chi_e*W_e + I + IH + (m_e*G_e*R_K*q_K(-1)*K(-1) + m_m*G_m*R_H*q_H(-1)*H_m(-1)+ mu_F*G_H*(Gamma_e - m_e*G_e)*R_K*q_K(-1)*K(-1) + mu_H*G_H*(Gamma_m - m_m*G_m)*R_H*q_H(-1)*H_m(-1)))/Y;
res_H = 100*(H - H_s - H_m)/H;


//************************************************************
// SHOCKS 
// Note that each of these shocks potentially has a component 
// that hits immediately and a 'news shock' component
//************************************************************
log(A)   = rhoA*log(A(-1)) - epsiA/100;
log(EJ)  = rhoJ*log(EJ(-1)) + epsiJ/100;
log(EK)  = rhoK*log(EK(-1)) + epsiK/100;
log(EH)  = rhoK*log(EH(-1)) + epsiH/100;
log(ESe) = rhoSe*log(ESe(-1)) + epsiSe/100;
log(ESm) = rhoSm*log(ESm(-1)) + epsiSm/100;
//Bank risk shocks are perfectly correlated across banks
log(ESF) = rhoSF*log(ESF(-1)) + epsiSF/100;
log(ESH) = rhoSH*log(ESH(-1)) + epsiSF/100;
log(EWe) = rhoWe*log(EWe(-1)) - epsiWe/100;
log(EWb) = rhoWb*log(EWb(-1)) - epsiWb/100;
//log(ERW) = rhoRW*log(ERW(-1))+epsiRW/100;
//ELTV = 1;
//Depreciation shocks are perfectly correlated across H and K sectors
EdH      = rhoHd*EdH(-1) + epsiHd/100;
EdK      = rhoHk*EdK(-1) + epsiHd/100;
end; 
//************************************************************
// Here specify the innovations
// Note that it is possible to shut down a shock by assigning it 
// a zero variance.
//************************************************************
//shocks;
//var epsiA   = sigma_epsiA;
//var epsiJ   = sigma_epsiJ;
//var epsiK   = sigma_epsiK;
//var epsiH   = sigma_epsiSH;
//var epsiSe  = sigma_epsiSe;
//var epsiSm  = sigma_epsiSm;
//var epsiSF  = sigma_epsiSF;
//var epsiSH  = sigma_epsiSH;
//var epsiWb  = sigma_epsiWb;
//var epsiWe  = sigma_epsiWe;
//var epsiHd  = sigma_epsiHd;
//var epsiHk  = sigma_epsiHk;
//var epsiRW  = sigma_epsiRW;
//end;

//*******************************
// SS
//*******************************
resid;
//steady;
//check;
//*******************************
// SIMULATIONS - IRFS
//*******************************

initval;
ELTV = 1;
//ERW=1;
end;

steady;
//resid;

endval;
ELTV = 1.05;
//ERW = 1.05;
end;

//steady;
//resid;
//simul(periods=1000);

yy = oo_.steady_state;
perfect_foresight_setup(periods=1000);
yy=[yy, oo_.endo_simul(:,2)];

perfect_foresight_solver;
yy=[yy, oo_.endo_simul(:,999)];
//ERW H_m_obs Y_net_obs C_obs I_obs IH_obs L_obs q_H_obs  q_K_obs b_m_obs  b_e_obs  def_rate_e  def_rate_m  av_def  RDsp  R_m_obs  R_F_obs  n_b_obs  n_e_obs  b_tot_obs W_m_obs



// Stochastic simulations (order=1 for linear approximation, order=2 for 2nd order approximation)
//stoch_simul(order=1,pruning,irf=40)
//Pr_H ERW H_m_obs Y_net_obs C_obs I_obs IH_obs L_obs q_H_obs  q_K_obs b_m_obs  b_e_obs  def_rate_e  def_rate_m  av_def  RDsp  R_m_obs  R_F_obs  n_b_obs  n_e_obs  b_tot_obs W_m_obs



//Vs Vm Ve Vb  def_rate_e def_rate_m def_rate_F def_rate_H x_e x_m
//Y R_D R_m R_H R_F H_m H_s b_m C C_m C_s D 
//E_F I K L L_m L_s n_b n_e q_H q_K r_K R_K
//R_tilde_F R_tilde_H rho_F rho_H Tr Tr_F Tr_H w
//Y_obs R_D_obs R_m_obs R_H R_F_obs H_m_obs H_s_obs b_m_obs C_obs C_m_obs C_s_obs D_obs 
//E_F_obs I_obs K_obs L_obs L_m_obs L_s_obs n_b_obs n_e_obs q_H_obs q_K_obs r_K_obs R_K_obs
//R_tilde_F_obs R_tilde_H_obs rho_F_obs rho_H_obs Tr_obs Tr_F_obs Tr_H_obs w_obs x_e x_m

//def_rate_e def_rate_m def_rate_F def_rate_H L_obs C_s_obs I_obs n_b_obs n_e_obs q_H_obs q_K_obs Y_obs Y_net_obs Y_net_2_obs C_m_obs R_m_obs R_F_obs R_D_obs x_e x_m
//R_tilde_H_obs R_tilde_F_obs R_K_obs b_m_obs b_e_obs H_m_obs H_s_obs bsp_H bsp_F phi_F phi_H deltaK deltaH res_chk rho_F_obs rho_H_obs EWe EWb

//def_rate_e def_rate_m def_rate_F def_rate_H L_obs C_s_obs I_obs n_b_obs n_e_obs q_H_obs q_K_obs Y_obs Y_net_obs C_m_obs R_m_obs R_F_obs R_D_obs x_e x_m
//R_tilde_H_obs R_tilde_F_obs b_m_obs b_e_obs H_m_obs H_s_obs bsp_H bsp_F phi_F phi_H rho_F_obs rho_H_obs 

//GDP_obs IH_obs H_obs L_s_obs L_m_obs 
//m_m m_e 
//av_def Y_net_obs Y_obs RDsp
//;   
//pruning,noprint, hp_filter=1600, noprint, hp_filter=1600, noprint,

//def_rate_e def_rate_m def_rate_F def_rate_H L_obs C_s_obs I_obs n_b_obs n_e_obs q_H_obs q_K_obs Y_obs Y_net_obs C_m_obs R_m_obs R_F_obs R_D_obs x_e x_m
//R_tilde_H_obs R_tilde_F_obs b_m_obs b_e_obs H_m_obs H_s_obs bsp_H bsp_F phi_F phi_H rho_F_obs rho_H_obs

