

var 
b_e      $B^e_t$ // entrepreneurial debt
g_H      $g^H_t$   // housing investment adjustment cost
g_H_1    // first derivative of g_H with respect to I_H
IH        $IH_t$  // housing investment 
PH       $PH_t$  // Profit of housing capital producing firm
A       $ A $   // productivity
b_m   $ b^m_t $     // mortgage debt
C       $ C_t $   // aggregate consumption
C_m     $ C^m_t $   // consumption of borrowers (impatient households)
C_s    $ C^s_t $    // consumption of savers (patient households)
D       $ D_t $   // aggregate deposits
def_rate_e $ Def^e_t $// corporate default rate
def_rate_m $ Def^m_t $// mortgage default rate//def_rate_B // default rate of corporate banks
def_rate_B $ Def^B_t $// default rate of mortgage banks//E_F        // equity invested in corporate banks

epsilonH $ \epsilon^H_t $
epsilonF $ \epsilon^F_t $
G_e    $ G^e_t $    // share of entrepreneurial capital belonging to firms that default (BGG parameter)
G_e_1      // First derivative of G_e with respect to omega_e (BGG parameter)//G_B        // share of corporate loans belonging to corporate banks that default (BGG parameter)
//G_B        // share of mortgage loans belonging to mortgage banks that default (BGG parameter)
F_pB    $ F^{pB}_t $    //cdf for bank default
G_B      $ G^B_t $  // share of mortgage loans belonging to mortgage banks that default (BGG parameter)
//F_pB
g_I     $ g^I_t $   // Capital investment adjustment cost 
g_I_1      // First derivative of g_I with respect to I_K
G_m     $ G^m_t $   // share of housing belonging to households that default (BGG parameter)
G_m_1      // First derivative of G_m with respect to omega_m (BGG parameter)
Gamma_e  $ \Gamma^e_t $  // Share of gross corporate revenues going to the bank (BGG parameter)
Gamma_e_1  // First derivative of Gamma_e with respect to omega_e(BGG parameter)//Gamma_B    // Share of gross corporate bank revenues going to depositors(BGG parameter)//Gamma_B_1  // First derivative of Gamma_B with respect to omega_F(BGG parameter)
Gamma_B $ \Gamma^b_t $
Gamma_B_1
//Gamma_B    // Share of gross mortgage bank revenues going to depositors(BGG parameter)
//Gamma_B_1  // First derivative of Gamma_B with respect to omega_H(BGG parameter)
Gamma_m $ \Gamma^m_t $   // Share of gross returns from housing going to the bank(BGG parameter)
Gamma_m_1  // First derivative of Gamma_m with respect to omega_m(BGG parameter)
F_pi    $ F^{pi}_t $   // cdf of probability of default
F_pe    $ F^{pe}_t $    //probability of defaultng households
//F_pi_1
//F_pe_1
H    $ H_t $      // Aggregate housing supply
H_m     $ H^m_t $   // Housing used by borrowers
H_s     $ H^s_t $   // Housing used by savers
I      $ I_t $    // Aggregate corporate investment
K     $ K_t $     // Aggregate capital stock 
L      $ L_t $    // Aggregate labour supply
L_m    $ L^m_t $    // Borrowers' labour supply
L_s     $ L^s_t $   // Savers' labour supply       
Lambda_m   $ \lambda^m_t $// LM on the budget constraint in the household borrower's problem
Lambda_s  $ \lambda^s_t $ // LM on the budget constraint in the household saver's problem
n_b    $ n^b_t $    // Bankers' net worth
n_e    $ n^e_t $    // Entrepreneurs' net worth
omega_bar_e $ \bar{\omega^e}_t $	// Idiosyncratic productivity shock below which the corporate borrower defaults
//omega_bar_B // Idiosyncratic corporate bank loan return shock below which the corporate bank defaults	
omega_bar_B	$ \bar{\omega^B}_t $// Idiosyncratic mortgage bank loan return shock below which the mortgage bank defaults
omega_bar_m $ \bar{\omega^m}_t $// Idiosyncratic housing return shock below which the mortgage borrower defaults	
PI      $ \PI     $   // Profit of business capital producing firm 
q_H      $ q^H_t  $    // Housing price
q_K     $  q^K_t $     // Capital price
R_D     $  R^D_t $     // Deposit interest rate
R_F    $ R^F_t  $      // Corporate loan interest rate NOMINAL
R_H    $   R^H_t$      // Aggregate financial return on housing (i.e. excluding imputed rents)
r_K      $ r^K_t  $    // Capital rental rate
R_K      $  R^K_t $    // Capital rate of return
R_m     $  R^m_t $     // Mortgage interest rate NOMINAL

R_tilde_F 	$  \tilde{R}^F_t $ // Aggregate return on a diversified corporate loan portfolio (i.e. portfolio return after accounting for loan losses)
R_tilde_H 	$  \tilde{R}^H_t $ // Aggregate return on a diversified housing loan portfolio (i.e. portfolio return after accounting for loan losses)
//rho_F       // Corporate bank return on equity
//rho_H       // Mortgage bank return on equity
//Tr_H
//Tr_F
Tr      $ Tr  $    // Total deposit insurance payments to banks//Tr_F        // deposit insurance payments to corporate banksTr_H        // deposit insurance payments to mortgage banks
UC_m    $ UC^m_t  $     // Borrowers' utility from consumption
UC_m_1      // Borrowers' marginal utility from consumption
UC_s    $ UC^s_t  $     // Savers' utility from consumption
UC_s_1      // Savers' marginal utility from consumption
UH_m  $ UC^h_t  $       // Borrowers' utility from housing services
UH_m_1      // Borrowers' marginal utility from housing services
UH_s    $  UH^s_t $     // Savers' utility from housing services
UH_s_1      // Savers' marginal utility from housing services
UL_m   $  UK^m_t $      // Borrowers' disutility from work
UL_m_1      // Borrowers' marginal disutility from work
UL_s   $  UL^s_t $      // Savers' disutility from work
UL_s_1      // Savers' marginal disutility from work
Util_m   $ Util^m_t  $    // Total borrower utility
Util_s   $  Util^s_t $    // Total saver utility
w     $ w_t  $       // Wage rate
W_b   $ W^b_t  $       // Bankers' wealth (pre-dividend)
Pr_H    $ Pr^H_t  $     //Profit of the bank
//rho         //return on bak equity
W_e    $ W^e_t  $      // Entrepreneurs' wealth (pre-dividend)
x_e    $ X^e_t  $      // Corporate leverage
x_m   $ X^m_t  $       // Household leverage
xi_e   $ \xi^e_t  $      // LM on bank's participation constraint in the entrepreneurs' problem
xi_m   $  \xi^m_t $      // LM on bank's participation constraint in the household borrower's problem
//xi_b          // LM on bank's BS constraint 
Y     $ Y_t  $       // Output
Y_obs    // All _obs variables are transformations to aid plotting. Transformations may differ across variables depending on what makes sense.
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
//R_tilde_F_obs
//R_tilde_H_obs
//Tr_obs
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
//ESB          // Shock (mortgage bank risk)
//ESB          // Shock (corporate bank risk)
ESB
EWe          // Shock (Entrepreneur net-worth)
EWb          // Shock (Banker net-worth)
EdH          // Shock (housing depreciation)
EdK          // Shock (capital depreciation)
markup_m
markup_F
EC
ECAB
EL
EbH
EbF
ELTVH
ELTVF
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
//rho_F
//rho_H
//rho_F_obs
//rho_H_obs
roe_B
DY_NET
DY_NET_OBS
credit_gap_obs
bsp_H_obs
bsp_F_obs
//variables for measurement equations
//credit_gap_obs
     //bsp_H_obs
       // bsp_F_obs
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
  //  bsp_H_data
 //  bsp_F_data
    dw_data
    dinve_data
    dc_data
    CR_data
welfare_
welfare_obs
dl_data
pb_ratio_data

;
parameters sigma_epsiHd sigma_epsiHk sigma_epsiA sigma_epsiJ sigma_epsiK sigma_epsiH 
sigma_epsiSe sigma_epsiSm sigma_epsiSB sigma_epsiWb sigma_epsiRW sigma_epsiWe sigma_epsiEC sigma_epsiEL pp  sigma_ELTVH sigma_ELTVF
hab Cyphi_H Cyphi_F phi_Fs phi_Hs alphaa delta_K delta_H betta_m betta_s mu_m mu_e mu_B
sigma_e1 sigma_m1 sigma_B varphi_s varphi_m v_s v_m chi_b chi_e eta  a_e a_s a_b   psi_i psi_h 
rhoA rhoJ rhoK rhoH rhoSe rhoSm rhoSB rhoWb rhoWe rhoHd rhoHk rhoRW rho_markup_m rho_markup_F rhoEC rhoECAB rhoEL rhoEbH rhoEbF rhoLTVH rhoLTVF
zeta_m zeta_F epsilonH1s epsilonF1s tau_H tau_F rp rpe phiinf1 kappa nu psib phis LTVHrule LTVFrule
gamma_y gamma_w gamma_inve gamma_c gamma_dbe gamma_dbm gamma_bspH gamma_dq_H def_rate_ss  omikronH omikronF pb_mean
;
    
varexo epsiA  epsiJ epsiK epsiSe epsiSm epsiSB epsiWb  epsiWe epsiH epsiHd epsiHk epsimarkup_m epsimarkup_F epsiEC epsiECAB epsiEL epsiEbH epsiEbF epsiLTVH epsiLTVF;

//===============================
//SET PARAMETER VALUES
load LTV1_parameter_values.mat;
//M=mortgage bank, E=entrepreneur, H=household, F=corporate bank
sigma_epsiHd=set_sigma_epsiHd;  //  shock to housing depreciation cost in EdH and EdK
sigma_epsiHk=set_sigma_epsiHk; // REDUNDANT, same as Hd in H and K sectorsf
 sigma_epsiA= set_sigma_epsiA; // productivity shock in A
 sigma_epsiJ= set_sigma_epsiJ; // housing preference shock in EJ
 sigma_epsiK = set_sigma_epsiK; // capital investment shock in EK
sigma_epsiH =set_sigma_epsiH ; // housing investment shock in EH
sigma_epsiSe=set_sigma_epsiSe; // entrepreneur risk shock in ESe
sigma_epsiSm=set_sigma_epsiSm; // housing risk shock in ESm
sigma_epsimarkup_m=set_sigma_epsimarkup_m;
sigma_epsimarkup_F=set_sigma_epsimarkup_F;
 // sigma_epsiSF=  set_sigma_epsiSF; // corporate bank risk shock in ESB===> why is this so large?R 
 //sigma_epsiSH= set_sigma_epsiSH; // mortgage bank risk shock is ESB
sigma_epsiSB=set_sigma_epsiSB;
  sigma_epsiWb=set_sigma_epsiWb;// banker net worth shock in EWb
  sigma_epsiRW=set_sigma_epsiRW;// REDUNDANT, RW shock in ERW? what is this? 
  sigma_epsiWe=  set_sigma_epsiWe; //  entrepreneur net worth shock in EWe
sigma_epsiEC=set_sigma_epsiEC;
sigma_epsiECAB=set_sigma_epsiECAB;
sigma_epsiEL=set_sigma_epsiEL;
sigma_epsiEbH=set_sigma_epsiEbH;
sigma_epsiEbF=set_sigma_epsiEbF;
sigma_ELTVH=set_sigma_ELTVH;
sigma_ELTVF=set_sigma_ELTVF;

 pp  = set_pp ; // transaction cost when recovering funds from failed banks
betta_s=set_betta_s; // patient HH discount factor
betta_m=set_betta_m; // impatient HH discount factor
hab =set_hab; // habit formation
Cyphi_H =set_Cyphi_H ; // REDUNDANT? parameter on capital requirement on mortgage banks
phi_Hs=set_phi_Hs; // REDUNDANT? same equation as above, equation shut off? 
Cyphi_F=set_Cyphi_F; // REDUNDANT? parameter on capital requirement on corporate banks
phi_Fs =set_phi_Fs  ; // REDUNDANT? appears in the same equation as above, the equation is shut off. 
alphaa=set_alphaa; // share of capital in output 
rp=set_rp; // loan repayment rate of impatient households
rpe=set_rpe; // loan repayment rate of entrepreneurs
delta_K= set_delta_K; // depreciation rate of housing & capital shocks
delta_H=set_delta_H; // same as above
mu_m=set_mu_m;
mu_e=set_mu_e;
mu_B=set_mu_B;
//mu_B=set_mu_B;
//mu_B=set_mu_B;
sigma_e1=set_sigma_e1;
sigma_m1=set_sigma_m1;
sigma_B=set_sigma_B;
//sigma_B=set_sigma_B;
//sigma_B=set_sigma_B;// all 8 are hyperparameters appearing in default cdfs?
varphi_s=set_varphi_s; //  patient HH preference parameter in utility
 varphi_m= set_varphi_m; // impatient HH preference parameter in utility
v_s=set_v_s; // patient HH preference parameter in utility 
v_m=set_v_m; // impatient HH preference parameter in utility
chi_b=set_chi_b; // banker preference parameter in utility
chi_e=set_chi_e; // entrepreneur preference parameter in utility 
  eta = set_eta; // patient HH  inverse frisch elasticity of labor supply 
 a_e =set_a_e; // ??? 
a_s=set_a_s; // ??? appears in the budget constraint of savers
 a_b = set_a_b; //  ??? appears in the budget constraint of borrowers
psi_i=set_psi_i; // parameter is capital adjustment cost function 
psi_h=set_psi_h; // same as above. Why are there two of these? 
rhoA =set_rhoA;
rhoJ =set_rhoJ ;
rhoK =set_rhoK;
 rhoSe = set_rhoSe;
rhoSm =set_rhoSm; 
//rhoSF =set_rhoSF; 
//rhoSH =set_rhoSH; 
rhoSB=set_rhoSB;
 rhoWb = set_rhoWb; 
rhoWe =set_rhoWe; 
rhoHd =set_rhoHd; 
rhoH=set_rhoH;
rhoHk =set_rhoHk ; 
rho_markup_m=set_rho_markup_m;
rho_markup_F=set_rho_markup_F;
rhoEC=set_rhoEC;
rhoECAB=set_rhoECAB;
rhoEL=set_rhoEL;
rhoEbH=set_rhoEbH;
rhoEbF=set_rhoEbF;
rhoLTVH=set_rhoLTVH;
rhoLTVF=set_rhoLTVF;

rhoRW =set_rhoRW ; //shock persistence parameters, same notation as standard deviation
zeta_m=set_zeta_m; // Interest rate stickiness, same as below, appears in  R_m R_F
zeta_F=set_zeta_F; // Interest rate stickiness,  appears in 51 FOC for bank business lending, 51 is not in use 
epsilonH1s=set_epsilonH1s;//this is the part of LTV rule that is not procyclical
epsilonF1s=set_epsilonF1s;// ??? =epsilonH  & =epsilonF respectively. What are these two? 
tau_H=set_tau_H; // appears in A1, B1, C1, D1, R_mi R_Fi
tau_F=set_tau_F;
phiinf1=set_phiinf1; //reaction to inflation--why is this so low? 
 kappa = set_kappa  ; // interest rate smoothing
nu=set_nu; // ??? appears in A1, C1
 psib = set_psib; //??? appears in A1, C1
phis=set_phis; /// ???? this is equated to phi, phi is not used anywhere??  is it the same as phi_F and phi_H maybe? 
LTVHrule=set_LTVHrule;
LTVFrule=set_LTVFrule;
gamma_y=set_gamma_y;
gamma_w=set_gamma_w;
gamma_inve=set_gamma_inve;
gamma_c=set_gamma_c;
gamma_dbe=set_gamma_dbe;
gamma_dbm=set_gamma_dbm;
gamma_bspH=set_gamma_bspH;
gamma_dq_H=set_gamma_dq_H;
def_rate_ss=set_def_rate_ss;
omikronH=set_omikronH;
omikronF=set_omikronF;
pb_mean=set_pb_mean;
//************************************************************
model;
//************************************************************
// Households 
//************************************************************
//*******************************
// Utility functions
//*******************************
//1
UC_s =EC* log(C_s-hab*C_s(-1)); 
//2
UC_m =EC*  log(C_m-hab*C_m(-1)); 
//3
UL_s = EL*varphi_s*L_s^(1+eta)/(1+eta); 
//4
UL_m = EL* varphi_m*L_m^(1+eta)/(1+eta); 
//5
UH_s = EJ*v_s*log(H_s(-1)); 
//6
UH_m = EJ*v_m*log(H_m(-1)); 
//7
Util_s = UC_s - UL_s + UH_s;
//8
Util_m = UC_m - UL_m + UH_m;
//9
UC_s_1 = EC* 1/(C_s-hab*C_s(-1));
//10
UC_m_1 = EC* 1/(C_m-hab*C_m(-1));
//11
UL_s_1 = EL*varphi_s*L_s^(eta);
//12
UL_m_1 = EL*varphi_m*L_m^(eta);
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
R_DD=R_D(-1)*(1-pp*(av_def));
// 23 foc H_s
Lambda_s*(q_H) = betta_s*UH_s_1(1) + betta_s*Lambda_s(1)*(1-deltaH(1))*EbH*q_H(1);
// 24 Budget Constraint
  C_s - C_e - C_b + q_H*(H_s-(1-deltaH)*H_s(-1))+ D = w*L_s +R_DD*D(-1) - Tr*a_s + PI + PH + (1-EWb)*(W_b)+(1-EWe)*(1-Gamma_e)*R_K*q_K(-1)*K(-1);
//*******************************
// Borrowers
//*******************************
//25 Default cut off borrowers depends on leverage
omega_bar_m = x_m(-1)/R_H;
//26 Household leverage definition
x_m = (R_m)*b_m/(H_m*q_H);
//27 Housing rate of return
R_H = (1-deltaH)*q_H/q_H(-1);
//28 foc C_m    
Lambda_m = UC_m_1;
//29 foc L_m
UL_m_1 = w*Lambda_m;
// 30--ltv limit borrowers
(b_m-b_m(-1)*(1-rp)*(1-F_pi))*R_m=ELTVH*(epsilonH-omikronH*phi_H)*(H_m-H_m(-1)*(1-deltaH))*EbH*q_H(+1);
// 31--ltv limit entrepreneurs
(b_e-b_e(-1)*(1-rpe)*(1-F_pe))*R_F=ELTVF*(epsilonF-omikronF*phi_F)*(K-K(-1)*(1-deltaK))*EbF*q_K(+1);
// 32--ltv rule borrowers
//epsilonH=epsilonH1s;
exp(epsilonH-epsilonH1s)=(b_m/steady_state(b_m))^(-LTVHrule);
//exp(epsilonH-epsilonH1s)=(bsp_H/bsp_H(-1))^(-LTVHrule);
// 33--ltv rule entrepreneurs
//epsilonF=epsilonF1s;
exp(epsilonF-epsilonF1s)=(b_e/steady_state(b_e))^(-LTVFrule);
//exp(epsilonF-epsilonF1s)=(bsp_F/bsp_F(-1))^(-LTVFrule);
// 34 foc H_m
betta_m*UH_m_1(1) - Lambda_m*(q_H) + betta_m*Lambda_m(1)*((1-G_m(1))*R_H(1)*q_H)+(xi_m)*(epsilonH*(EbH*q_H(+1)))-xi_m(+1)*(epsilonH*(EbH(+1)*q_H(+2)))*(1-deltaH) = 0;  //needs to be changed                                       
// 35 foc b_m 
Lambda_m-betta_m*Lambda_m(1)*(R_m*(1-F_pi(+1)))-xi_m*R_m+xi_m(+1)*(1-rp)*(1-F_pi(+1))*R_m(+1)=0;///needs to be changed
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
x_e = (R_F)*(q_K*K-n_e)/(q_K*K);
//  FOC x_e
//- Gamma_e_1(1) + xi_e*((1-Gamma_B(1))*(Gamma_e_1(1)- m_e*G_e_1(1)));             
// 40 Foc k
(((1-G_e(1))*R_K(1)*q_K))*betta_s*(Lambda_s(+1)/Lambda_s)-(q_K*R_F*(1-F_pe(+1)))+xi_e*(epsilonF*EbF*q_K(+1)-R_F)-xi_e(+1)*(epsilonF*EbF(1)*q_K(+2)*(1-deltaK)-(1-rpe)*(1-F_pe(+1))*R_F(+1))=0 ;//need to be changed for pc
// 41 Wealth dynamics entrepreneurs  
W_e= EWe*(1-Gamma_e)*R_K*q_K(-1)*K(-1) ;
// 42 Evolution individual net worth 
n_e = (1-chi_e)*W_e;
// 43 Corporate dividends
C_e=chi_e*W_e;
//***********************************************************
// Bankers
//************************************************************
// 42 banker wealth
W_b=EWb*Pr_H;
// 45 Banker net worth allocated to lending
n_b = (1-chi_b)*W_b;
// 46 Bank dividends
C_b=chi_b*(W_b);
//************************************************************
// Banks
//************************************************************
// Definitions
// 47 Default threshold bank
//omega_bar_B=(1-phib)*R_D(-1)/((b_m*R_tilde_H+b_e*R_tilde_F)/(b_e+b_m));
omega_bar_B=R_D(-1)*D(-1)/(R_tilde_H*b_m+R_tilde_F*b_e);
// 49 bank profit
Pr_H(+1) =
(1-G_B(+1))*((1-F_pi(+1))*b_m*R_m+G_m(+1)*(1 - mu_m)*(b_m*(R_m)/omega_bar_m(+1)))
+(1-G_B(+1))*(((1-F_pe(+1))*b_e*R_F+G_e(+1)*(1 - mu_e)*(b_e*(R_F)/omega_bar_e(+1))))
-(1-F_pB(+1))*(b_m*(1-phib))*R_D
-(1-F_pB(+1))*(b_e*(1-phib))*R_D;
//50 FOC for bank mortgage lending
//capital ratio
phib=ECAB*n_b/((((R_mi/R_m)^(-tau_H))*b_m+((R_Fi/R_F)^(-tau_F))*b_e));
//phib=ECAB*n_b/(b_m+b_e);
//phib=n_b/(b_m+b_e);
//FOC wrt R_Fi
A1 = (betta_s*zeta_m)*(Lambda_s(+1)*((1-F_pB(+1))*(R_D)+nu*(phib/phi_H)^(1-psib)/(((R_mi/R_m)^(-tau_H))*b_m+((R_Fi/R_F)^(-tau_F))*b_e))*((R_m)^(tau_H))*b_m)+(betta_s*zeta_m)*A1(+1);
B1=(betta_s*zeta_m)*Lambda_s(+1)*((1-G_B(+1))*(((1-F_pi(+1)))+(G_m(+1)*(1 - mu_m)/omega_bar_m(+1))))*((((R_m)^(tau_H))*b_m))+(betta_s*zeta_m)*B1(+1);
//FOC wrt R_mi
C1=(betta_s*zeta_F)*(Lambda_s(+1)*((1-F_pB(+1))*(R_D)+(nu*(phib/phi_F)^(1-psib)/(((R_mi/R_m)^(-tau_H))*b_m+((R_Fi/R_F)^(-tau_F))*b_e)))*((R_F)^(tau_F))*b_e)+(betta_s*zeta_F)*C1(+1);
D1=(betta_s*zeta_F)*Lambda_s(+1)*((1-G_B(+1))*(((1-F_pe(+1)))+(G_e(+1)*(1 - mu_e)/omega_bar_e(+1))))*((((R_F)^(tau_F))*b_e))+(betta_s*zeta_F)*D1(+1);
R_mi=markup_m*(tau_H/(tau_H-1))*(A1/(B1));
R_Fi=markup_F*(tau_F/(tau_F-1))*(C1/(D1));
R_m=((1-zeta_m)*R_mi^(1-tau_H)+zeta_m*R_m(-1)^(1-tau_H))^(1/(1-tau_H));
R_F=((1-zeta_F)*R_Fi^(1-tau_F)+zeta_F*R_F(-1)^(1-tau_F))^(1/(1-tau_F));
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
g_I = psi_i/2*(I/I(-1)-1)^2;
// 60 foc of (57) wrt I
g_I_1 = psi_i*(I/I(-1)-1);
// 61 Capital stock evolution
K = (1-deltaK)*K(-1) + I*(1-g_I );
// 62 Capital price
q_K/EK = 1 + g_I + I/I(-1)*g_I_1 - betta_s*Lambda_s(1)/Lambda_s*(I(1)/I)^2*g_I_1(1);
// 63 Flow profit of capital producers
PI = q_K*I - (1 + g_I)*I; 
//*******************************
// Inside equity market
//*******************************
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
g_H = psi_h/2*(IH/IH(-1)-1)^2;
// 68 foc of (63) wrt IH
g_H_1 = psi_h*(IH/IH(-1)-1);
// 69 Housing stock evolution
H = (1-deltaH)*H(-1) + IH*(1-g_H);
// 70 Housing price with CEE Adj costs
q_H/EH = 1 + g_H + IH/IH(-1)*g_H_1 - betta_s*Lambda_s(1)/Lambda_s*(IH(1)/IH)^2*g_H_1(1);
// 71 Flow profits by housing producers 
PH= q_H*IH - (1 + g_H)*IH;  
// 72 Aggregate housing market clearing (Supply (LHS) = Demand (RHS))
H = H_m + H_s; 
//*******************************
// Deposit insurance agency
//*******************************
//  Deposit insurance transfers to corporate bank depositors whose banks have gone bankrupt
//Tr_F = (omega_bar_B - Gamma_B + mu_B*G_B)*R_tilde_F*(q_K(-1)*K(-1)-(1-chi_e)*W_e(-1));
// 70 Deposit insurance transfers to household bank depositors whose banks have gone bankrupt
//Tr_H = (omega_bar_B - Gamma_B + mu_B*G_B)*((((1-F_pi)*b_m(-1)*R_m(-1)+G_m*(1 - mu_m)*(b_m(-1)*(R_m(-1))/omega_bar_m))));
//Tr_F = (omega_bar_B - Gamma_B + mu_B*G_B)*(((1-F_pe)*b_e(-1)*R_F(-1)+G_e*(1 - mu_e)*(b_e(-1)*(R_F(-1))/omega_bar_e)));
// 73
// 71 Total deposit insurance transfers
//Tr = Tr_F + Tr_H;
Tr=(omega_bar_B - Gamma_B + mu_B*G_B)*((((1-F_pi)*b_m(-1)*R_m(-1)+(1-F_pe)*b_e(-1)*R_F(-1)+
G_m*(1 - mu_m)*(b_m(-1)*(R_m(-1))/omega_bar_m)+G_e*(1 - mu_e)*(b_e(-1)*(R_F(-1))/omega_bar_e))));
//**************************************************************
//Taylor rule
//**************************************************************
//************************************************************
// Distribution functions
//************************************************************
// 74
Gamma_m   = normcdf((log(omega_bar_m)- ((ESm)*sigma_m1)^2/2)/((ESm)*sigma_m1)) + omega_bar_m*(1-normcdf((log(omega_bar_m)+(ESm)*sigma_m1^2/2)/((ESm)*sigma_m1)));
// 75
Gamma_e   = normcdf((log(omega_bar_e) - ((ESe)*sigma_e1)^2/2)/((ESe)*sigma_e1)) + omega_bar_e*(1-normcdf((log(omega_bar_e)+(ESe)*sigma_e1 ^2/2)/((ESe)*sigma_e1)));
// 76
//Gamma_B   = normcdf((log(omega_bar_B) - ((ESB)*sigma_B)^2/2)/((ESB)*sigma_B)) + omega_bar_B*(1-normcdf((log(omega_bar_B)+(ESB)*sigma_B^2/2)/((ESB)*sigma_B)));
// 77
Gamma_B   = normcdf((log(omega_bar_B)-((ESB)*sigma_B)^2/2)/((ESB)*sigma_B)) + omega_bar_B*(1-normcdf((log(omega_bar_B)+(ESB)*sigma_B^2/2)/((ESB)*sigma_B)));
// 78
Gamma_m_1 = (normpdf((log(omega_bar_m)-(ESm*sigma_m1)^2/2)/(ESm*sigma_m1))-omega_bar_m*normpdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1)))/(ESm*sigma_m1*omega_bar_m) + (1-normcdf((log(omega_bar_m)+(ESm*sigma_m1)^2/2)/(ESm*sigma_m1)));
// 79
Gamma_e_1 = (normpdf((log(omega_bar_e)-(ESe*sigma_e1) ^2/2)/(ESe*sigma_e1))-omega_bar_e*normpdf((log(omega_bar_e)+ESe*sigma_e1 ^2/2)/(ESe*sigma_e1)))/(ESe*sigma_e1 *omega_bar_e) + (1-normcdf((log(omega_bar_e)+(ESe*sigma_e1)^2/2)/(ESe*sigma_e1)));
// 80
//Gamma_B_1 = (normpdf((log(omega_bar_B)-(ESB*sigma_B)^2/2)/(ESB*sigma_B))-omega_bar_B*normpdf((log(omega_bar_B)+ESB*sigma_B^2/2)/(ESB*sigma_B)))/(ESB*sigma_B*omega_bar_B) + (1-normcdf((log(omega_bar_B)+(ESB*sigma_B)^2/2)/(ESB*sigma_B)));
// 81
Gamma_B_1 = (normpdf((log(omega_bar_B)-(ESB*sigma_B)^2/2)/(ESB*sigma_B))-omega_bar_B*normpdf((log(omega_bar_B)+ESB*sigma_B^2/2)/(ESB*sigma_B)))/(ESB*sigma_B*omega_bar_B) + (1-normcdf((log(omega_bar_B)+(ESB*sigma_B)^2/2)/(ESB*sigma_B)));
// 82
G_m   = normcdf((log(omega_bar_m)-ESm*sigma_m1^2/2)/(ESm*sigma_m1));
// 83
G_e   = normcdf((log(omega_bar_e)-ESe*sigma_e1 ^2/2)/(ESe*sigma_e1));
// 84
//G_B   = normcdf((log(omega_bar_B)-ESB*sigma_B^2/2)/(ESB*sigma_B));
// 85
G_B   = normcdf((log(omega_bar_B)-ESB*sigma_B^2/2)/(ESB*sigma_B));
// 86
G_m_1 = normpdf((log(omega_bar_m)-ESm*sigma_m1^2/2)/(ESm*sigma_m1))/(ESm*sigma_m1*omega_bar_m);
// 87
F_pi = normcdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1));//probability of defaultng households
//F_pi_1 = normpdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1))/(ESm*sigma_m1*omega_bar_m);
// 95
F_pe = normcdf((log(omega_bar_e)+ESe*sigma_e1^2/2)/(ESe*sigma_e1));//probability of defaultng households
//F_pe_1 = normpdf((log(omega_bar_e)+ESm*sigma_e1^2/2)/(ESm*sigma_e1))/(ESe*sigma_e1 *omega_bar_e);//probability of defaultng households
// 96
F_pB   = normcdf((log(omega_bar_B)+ESB*sigma_B^2/2)/(ESB*sigma_B));
// 97
//F_pB   = normcdf((log(omega_bar_B)+ESB*sigma_B^2/2)/(ESB*sigma_B));
// 91
G_e_1 = normpdf((log(omega_bar_e)-ESe*sigma_e1 ^2/2)/(ESe*sigma_e1))/(ESe*sigma_e1 *omega_bar_e);
// 92 Capital requirement on mortgage banks
//exp(phi_H-phi_Hs) = (b_m/b_m(-1))^(Cyphi_H);
exp(phi_H-phi_Hs) = (b_m/steady_state(b_m))^(Cyphi_H);//original version
//exp(phi_H-phi_Hs)=((Y_net_2)/(Y_net_2(-1)))^(Cyphi_H);
// 93 Capital requirement on corporate banks
//exp(phi_F-phi_Fs) = (b_e/b_e(-1))^(Cyphi_F);
exp(phi_F-phi_Fs) = (b_e/steady_state(b_e))^(Cyphi_F);//original version
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
Y_net_2 = Y - (R_D(-1)*pp*av_def*D(-1)
+m_e*G_e*R_K*q_K(-1)*K(-1) 
+ m_m*G_m*R_H*q_H(-1)*H_m(-1)+
 mu_B*G_B*(Gamma_e - m_e*G_e)*R_K*q_K(-1)*K(-1) +
 mu_B*G_B*(Gamma_m - m_m*G_m)*R_H*q_H(-1)*H_m(-1));
GDP = Y + UL_s_1*H_s/Lambda_s + UL_m_1*H_m/Lambda_m;
// Monitoring cost
m_e = mu_e;
m_m = mu_m;
//Average default of banks
//av_def = ((1-phi_F)*b_e*def_rate_B + (1-phi_H)*b_m*def_rate_B)/D;
//av_def = ((1-phi_F)*b_e*def_rate_B + (1-phi_H)*b_m*def_rate_B)/((1-phi_F)*b_e+(1-phi_H)*b_m);
av_def=def_rate_B;
// Deposit spread
RDsp = 400*(R_D-R_DD(+1));
// Borrowers wealth
W_m = (1-Gamma_m)*R_H*q_H(-1)*H_m(-1) + w*L_m;
//************************************************************
// Reporting adjunct variables
//************************************************************
//def_rate_m = normcdf((log(omega_bar_m)+(ESm)*sigma_m1^2/2)/((ESm)*sigma_m1),0,1)*400;
//def_rate_e = normcdf((log(omega_bar_e)+(ESe)*(sigma_e1)^2/2)/((ESe)*(sigma_e1)),0,1)*400;
//def_rate_B = normcdf((log()+f(ESB)*sigma_B^2/2)/((ESB)*sigma_B))*400;
//def_rate_B = normcdf((log(omega_bar_B)+(ESB)*sigma_B^2/2)/((ESB)*sigma_B))*400;
def_rate_m = def_rate_ss+normcdf((log(omega_bar_m)+((ESm))*sigma_m1^2/2)/(((ESm))*sigma_m1))*400;
def_rate_e = def_rate_ss+normcdf((log(omega_bar_e)+((ESe))*(sigma_e1)^2/2)/(((ESe))*(sigma_e1)))*400;
//def_rate_B = def_rate_ss+normcdf((log(omega_bar_B)+((ESB))*sigma_B^2/2)/(((ESB))*sigma_B))*400;
def_rate_B = def_rate_ss+normcdf((log(omega_bar_B)+((ESB))*sigma_B^2/2)/(((ESB))*sigma_B))*400;
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
//R_tilde_F_obs = (R_tilde_F-steady_state(R_tilde_F))*400;
//R_tilde_H_obs = (R_tilde_H-steady_state(R_tilde_H))*400;
//Tr_obs = log(Tr/steady_state(Tr))*100;
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
res_chk = (R_D(-1)*pp*av_def*D(-1)+C + chi_b*W_b + chi_e*W_e + I + IH + (m_e*G_e*R_K*q_K(-1)*K(-1) + m_m*G_m*R_H*q_H(-1)*H_m(-1)+ mu_B*G_B*(Gamma_e - m_e*G_e)*R_K*q_K(-1)*K(-1) + mu_B*G_B*(Gamma_m - m_m*G_m)*R_H*q_H(-1)*H_m(-1)))/Y;
res_H = 100*(H - H_s - H_m)/H;
//************************************************************
// SHOCKS 
// Note that each of these shocks potentially has a component 
// that hits immediately and a 'news shock' component
//************************************************************
log(A)   = rhoA*log(A(-1))+ epsiA;
log(EJ)  = rhoJ*log(EJ(-1)) + epsiJ;
log(EK)  = rhoK*log(EK(-1)) + epsiK;
log(EH)  = rhoH*log(EH(-1)) + epsiH;
log(ESe) = rhoSe*log(ESe(-1)) + epsiSe ;
log(ESm) = rhoSm*log(ESm(-1)) + epsiSm;
//Bank risk shocks are perfectly correlated across banks--> this also has been changed it seems?
//log(ESB) = rhoSF*log(ESB(-1)) + epsiSF;
//log(ESB) = rhoSH*log(ESB(-1)) + epsiSH;
log(ESB) = rhoSB*log(ESB(-1)) + epsiSB;
log(markup_m)=rho_markup_m*log(markup_m(-1))+epsimarkup_m;
log(markup_F)=rho_markup_F*log(markup_F(-1))+epsimarkup_F;
log(EC)=rhoEC*log(EC(-1))+epsiEC;
log(ECAB)=rhoEC*log(ECAB(-1))+epsiECAB;
log(EL)=rhoEL*log(EL(-1))+epsiEL;
log(EbH)=rhoEbH*log(EbH(-1))+epsiEbH;
log(EbF)=rhoEbF*log(EbF(-1))+epsiEbF ;
log(EWe) = rhoWe*log(EWe(-1)) + epsiWe;
log(EWb) = rhoWb*log(EWb(-1)) + epsiWb ;
//log(ERW) = rhoRW*log(ERW(-1))+ epsiRW;
log(ELTVH) = rhoLTVH*log(ELTVH(-1))+epsiLTVH;
log(ELTVF) = rhoLTVF*log(ELTVF(-1))+epsiLTVF;
//Depreciation shocks are perfectly correlated across H and K sectors--> this has been changed?
EdH      = rhoHd*EdH(-1) + epsiHd;
EdK      = rhoHk*EdK(-1) + epsiHk;
// 49 Rate of return on corporate bank equity
//rho_F = (1-Gamma_B)*R_tilde_F/phi_F(-1);
// 50 Rate of return on mortgage bank equity
//rho_H =(1-Gamma_B)*R_tilde_H/phi_H(-1);
//rho_F_obs = (rho_F-steady_state(rho_F))*400;
//rho_H_obs = (rho_H-steady_state(rho_H))*400;
roe_B=((1-Gamma_B)*(R_tilde_F*b_e+R_tilde_H*b_m)/(b_e+b_m))/phib;
//DY_NET=0;
DY_NET_OBS=0;
DY_NET=Y_net-Y_net(-1);
//DY_NET_OBS=log(DY_NET/steady_state(DY_NET))*100;
welfare_=((C_m)*Vm+(C_s)*Vs)/((C_m)+(C_s));
welfare_obs=100*log(welfare_/steady_state(welfare_));
credit_gap_obs=100*log( (b_tot/steady_state(b_tot)) / (Y_net/steady_state(Y_net)));
//bsp_H_obs=100*log(bsp_H/steady_state(bsp_H));
bsp_H_obs=400*log((R_tilde_H(1)-R_D)/(steady_state(R_tilde_H)-steady_state(R_D)));
bsp_F_obs=400*log((R_tilde_F(1)-R_D)/(steady_state(R_tilde_F)-steady_state(R_D)));
//MEASUREMENT EQUATIONS
dy_data=Y_net_obs-Y_net_obs(-1)+gamma_y ;
//pb_ratio_data-pb_ratio_data(-1) = n_b_obs - n_b_obs(-1) ;
pb_ratio_data-pb_ratio_data(-1)= (1/ESB-1);
//dy_data=Y_net_2_obs-Y_net_2_obs(-1)+gamma_y;
dq_H_data=(q_H_obs-q_H_obs(-1))+gamma_dq_H ;
dbe_data=4*(b_e_obs-b_e_obs(-1))+gamma_dbe;
//dbe_data=100*log(b_e/b_e(-1));F
dbm_data=4*(b_m_obs-b_m_obs(-1))+gamma_dbm;
db_total_data=5.73 + 4*(b_m*(b_m_obs -b_m_obs(-1) ) + b_e*(b_e_obs - b_e_obs(-1)) )/(b_e+b_m);
bank_rate_data=R_D_obs+100*(steady_state(R_D)-1);
b_to_Y_data-b_to_Y_data(-1)=credit_gap_obs-credit_gap_obs(-1);
d_b_to_Y_data=b_to_Y_data-b_to_Y_data(-1);
int_rate_HH_data=R_m_obs+100*(steady_state(R_m)-1) ;
int_rate_business_data=R_F_obs+100*(steady_state(R_F)-1);
//d_bsp_H_data=(bsp_H_obs-bsp_H_obs(-1))/4;
//bsp_H_data-bsp_H_data(-1)=log((R_tilde_H(1)-R_D)/(R_tilde_H-R_D(-1)))+gamma_bspH;
//bsp_F_data-bsp_F_data(-1)=log((R_tilde_F(1)-R_D)/(R_tilde_F-R_D(-1)))+gamma_bspH;
dw_data=w_obs-w_obs(-1)+gamma_w;
//dinve_data=(I_obs-I_obs(-1))+(IH_obs-IH_obs(-1))+gamma_y;
dinve_data=I_obs-I_obs(-1)+(IH_obs-IH_obs(-1))+ gamma_inve ;
dc_data=C_obs-C_obs(-1)+gamma_c ;
CR_data=100*phi_H;
dl_data=0.095 + L - L(-1);
end; 
//************************************************************
// Here specify the innovations
// Note that it is possible to shut down a shock by assigning it 
// a zero variance.
//************************************************************
shocks;
var epsiA   = sigma_epsiA;
var epsiJ   = sigma_epsiJ;
var epsiK   = sigma_epsiK;
var epsiH   = sigma_epsiH;
var epsiSe  = sigma_epsiSe;
var epsiSm  = sigma_epsiSm;
//var epsiSF  = sigma_epsiSF;
//var epsiSH  = sigma_epsiSH;
var epsiSB = sigma_epsiSB;
var epsiWb  = sigma_epsiWb;
var epsiWe  = sigma_epsiWe; 
var epsimarkup_m= sigma_epsimarkup_m;
var epsimarkup_F=sigma_epsimarkup_F;
var epsiEC = sigma_epsiEC;
var epsiECAB=sigma_epsiECAB;
var epsiEL=sigma_epsiEL;
var epsiEbH=sigma_epsiEbH;
var epsiEbF=sigma_epsiEbF;
var epsiHd  = sigma_epsiHd;
var epsiHk  = sigma_epsiHk;
//var epsiRW  = sigma_epsiRW;
end;
//resid;
//steady;
//check;
estimated_params;
//CALIBRATION BLOCK
//delta_H, 0.01 ,0,1,BETA_PDF,0.05,0.025;
//delta_K,0.03 ,0,1,BETA_PDF,0.05,0.025;
//gamma_y,  0.21 ,1e-6,5,UNIFORM_PDF,0.5,  0.288675134594813;
//gamma_w,  0.8 ,1e-6,5,UNIFORM_PDF,0.5,  0.288675134594813;
//gamma_inve,  0.17 ,1e-6,5,UNIFORM_PDF,0.5,  0.288675134594813;
//gamma_c,  0.24 ,1e-6,5,UNIFORM_PDF,0.5,  0.288675134594813;
//chi_e, 0.0760,0,1,UNIFORM_PDF,0.5, 0.288675134594813;
//chi_b, 0.005 ,0.01,1,UNIFORM_PDF,0.5, 0.288675134594813;
//alphaa,0.3 ,0,1,BETA_PDF,0.3,0.05;
//mu_m, 0.3 ,0,1,BETA_PDF,0.3,0.05;
//mu_e, 0.3,0,1,BETA_PDF,0.3,0.05;
//mu_B,0.3,0,1,BETA_PDF,0.3,0.05;
//mu_B, 0.3,0,1,BETA_PDF,0.3,0.05;
//v_s,0.01  ,0.001,1,BETA_PDF,0.5,0.25;
//v_m,  0.076 ,0.01,1,BETA_PDF,0.5,0.25;
//pp,0.1,0,1,BETA_PDF,0.1,0.05;
//varphi_m, 0.7290 ,0,5,GAMMA_PDF,1,0.5;
//varphi_s,   1.049981184116970  ,0,5,GAMMA_PDF,1,0.5;
//tau_H,100,0,200,UNIFORM_PDF,100,      28.867513459481287;
//tau_F,100,0,200,UNIFORM_PDF,100,      28.867513459481287;
//ESTIMATION BLOCK
/*
rp,0.0100,0,1,BETA_PDF,0.05,0.025;
rpe,     0.0456,0,1,BETA_PDF,0.05,0.025;//danmarks national bank guy has values ~0.25 for these parameters.
betta_m,   0.9719 ,0.95,0.995,BETA_PDF,0.98,0.01;
hab,   0.6626,0,1,BETA_PDF,0.5,0.2;
*/
//sigma_e1,0.0459 ,1e-5,1,   INV_GAMMA_PDF,0.1,2;  //    UNIFORM_PDF,0.5,2.886751345948129;
//sigma_m1, 0.0335,1e-5,1,   INV_GAMMA_PDF,0.1,2;   //   UNIFORM_PDF,0.5,2.886751345948129;// GAMMA_PDF,0.3,0.05; //
//sigma_B,0.0658,1e-5,1,     INV_GAMMA_PDF,0.1,2;  //  UNIFORM_PDF,0.5,2.886751345948129; //  GAMMA_PDF,0.04,0.01;// 
psi_i, 8.4488  ,1,100,NORMAL_PDF,5,2;//standard prior from SW
psi_h,  10.2520 ,1,100,NORMAL_PDF,5,2;//standard prior from SW
//zeta_m,  0.7541  ,0,1,BETA_PDF,0.5,0.2;//standard prior for stickiness terms-->check whether another paper attempts to estimate this
//zeta_F,0.5042,0,1,BETA_PDF,0.5,0.2;//standard prior for stickiness terms-->check whether another paper attempts to estimate this
//stderr epsiHd,0.1,0,99,          UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
//stderr epsiHk,0.001,0,99,    UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
stderr epsiA,0.1,0,99,     UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
stderr epsiJ,0.1,0,99,   UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
//stderr epsiK,0.1,0,99,  UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
stderr epsiH,0.1,0,99,    UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
stderr epsiSe,0.1,0,99,  INV_GAMMA_PDF,0.1,2;//    UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
//stderr epsiSm,0.1,0,99,  INV_GAMMA_PDF,0.1,2;//    UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
stderr epsiSB,0.1,0,99,  INV_GAMMA_PDF,0.1,2;//    UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
stderr epsiWe,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
//stderr epsiWb,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
stderr epsimarkup_m,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
stderr epsimarkup_F,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
stderr epsiEC,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
stderr epsiECAB,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
stderr epsiLTVH,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
//stderr epsiLTVF,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//



//stderr epsiEL,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
//stderr epsiEbH,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
stderr epsiEbF,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
//stderr meas_dw,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
//stderr meas_dinve,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
//stderr meas_dc,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
//stderr epsiFirm,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
//stderr epsiBank,0.1,0,99, UNIFORM_PDF,10, 5.773502691896258;//  GAMMA_PDF,4,4;//
//rhoHd,0.01,0,1,BETA_PDF,0.5,0.2;
//rhoHk,0.01,0,1,BETA_PDF,0.5,0.2;
rhoA,0.01,0,1,BETA_PDF,0.5,0.2;
rhoJ,0.01,0,1,BETA_PDF,0.5,0.2;
//rhoK,0.01,0,1,BETA_PDF,0.5,0.2;
rhoH,0.01,0,1,BETA_PDF,0.5,0.2;
rhoSe,0.01,0,1,BETA_PDF,0.5,0.2;
//rhoSm,0.01,0,1,BETA_PDF,0.5,0.2;                     
rhoSB,0.01,0,1,BETA_PDF,0.5,0.2;                                                                                                                                       
rhoWe,0.01,0,1,BETA_PDF,0.5,0.2;
//rhoWb,0.01,0,1,BETA_PDF,0.5,0.2;
rho_markup_m,0.01,0,1,BETA_PDF,0.5,0.2;
rho_markup_F,0.01,0,1,BETA_PDF,0.5,0.2;
rhoEC,0.01,0,1,BETA_PDF,0.5,0.2;
rhoECAB,0.01,0,1,BETA_PDF,0.5,0.2;
//rhoEL,0.01,0,1,BETA_PDF,0.5,0.2;
//rhoEbH,0.01,0,1,BETA_PDF,0.5,0.2;
rhoEbF,0.01,0,1,BETA_PDF,0.5,0.2;
rhoLTVH,0.01,0,1,BETA_PDF,0.5,0.2;
//rhoLTVF,0.01,0,1,BETA_PDF,0.5,0.2;
end;
varobs
dy_data,dw_data,dc_data,dinve_data,
int_rate_HH_data,int_rate_business_data,bank_rate_data,dq_H_data,dbm_data,dbe_data,pb_ratio_data;//,
//,,;//,dbm_data,dbe_data;dl_data;,dl_data
//-->dc_data should be taken out potentially    ,dinve_data,dw_data,;,bsp_H_data,,dinve_data,;,dc_data, dbm_data;,, ;int_rate_HH_data , bank_rate_data ,dbe_data bank_rate_data  , int_rate_HH_data;bank_rate_data  ,dbe_data; dbe_data,dbm_data bsp_H_data;; , //////////////dq_H_data, b_to_Y_data, int_rate_HH_data, bank_rate_data,bsp_H_data;,b_to_Y_data ,b_to_Y_data,bsp_H_data
//bsp_H_data,
estimation(datafile='estimation_dataset_quarterly.mat', 
mode_compute=0,
//mode_check,
mode_file='LTV1_mode.mat',
//optim=('TolFun',1e-4,'TolX',1e-4,'nclimb_mh',20000,'ncov_mh',2000,'NumberOfMh',1),//,'UseParallel',1
//optim=('TolFun',1e-4),
optim=('TolFun',1e-7,'TolX',1e-7,'Display','iter','MaxIter',500,'Hessian','bfgs'),//,'UseParallel',
cova_compute=0,
//optim=('TolFun',1e-4,'MaxIter',20,'Hessian',2),//,'UseParallel',1
//optimset=('TolFun',1e-5,'TolX',1e-5,'Display','iter','MaxIter',50),
//optim=('TolFun',10e-4,'TolX',10e-4,'Display','iter','UseParallel',1,'MaxIter',9999999,'MaxFunEvals',9999),
kalman_algo=4,
diffuse_kalman_tol=1e-15,
//use_univariate_filters_if_singularity_is_detected=1,
//lik_init=1,
nograph,nodiagnostic,
first_obs=3,
nobs=72,
presample=8,
mh_nblocks=1,
//mcmc_jumping_covariance=prior_variance,
//load_mh_file,
mh_replic=0000,
//bayesian_irf,
//load_results_after_load_mh,
mh_jscale=0.2,
mh_drop=0.5);
//mh_init_scale=0.001);


//varexo epsiA  epsiJ epsiK epsiSe epsiSm epsiSF epsiSH epsiWb  epsiWe epsiH epsiHd epsiHk;


shock_groups;
'Productivity'=epsiA;
'Housing Preference'=epsiJ;
'House Price Shock'=epsiH;
'Business Risk'=epsiSe;
'Bank Risk'=epsiSB;
'Entrepreneur Net Worth'=epsiWe;
//'Banker Net Worth'=epsiWb;
'Mortgage Lending Mark-up'=epsimarkup_m;
'Corporate Lending Mark-up'=epsimarkup_F;
'Consumption Preference'=epsiEC;
'Bank Capital'=epsiECAB;
'Expected Capital Price'=epsiEbF;
//'Expected House Price'=epsiEbH;
'Household LTV'=epsiLTVH;
end;

 //stoch_simul(order=1,irf=100,nograph,nocorr,nomoments,nodecomposition,nodisplay,noprint,tex,relative_irf);
stoch_simul(order=1,irf=40);
shock_decomposition(first_obs=1,datafile='estimation_dataset_quarterly.mat',parameter_set=posterior_mode,use_shock_groups,init_state=0)
dy_data,dq_H_data,int_rate_HH_data,int_rate_business_data,bank_rate_data,dbe_data,dbm_data,dinve_data,dw_data,dc_data,pb_ratio_data;
