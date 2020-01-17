function set_parameter_values(data)
set_sigma_epsiHd= data.parameters.sigma_epsiHd ;  %  shock to housing depreciation cost in EdH and EdK %1
set_sigma_epsiHk=data.parameters.sigma_epsiHk; % shock to capital depreciation 

 set_sigma_epsiA= data.parameters.sigma_epsiA; % productivity shock in A %3
 set_sigma_epsiJ=data.parameters.sigma_epsiJ ; % housing preference shock in EJ %4
 
 set_sigma_epsiK =   data.parameters.sigma_epsiK; % capital investment shock in EK %5
set_sigma_epsiH = data.parameters.sigma_epsiH; % housing investment shock in EH %6

set_sigma_epsiSe=data.parameters.sigma_epsiSe; % entrepreneur risk shock in ESe %7
set_sigma_epsiSm=data.parameters.sigma_epsiSm; % housing risk shock in ESm %8fd
  %set_sigma_epsiSF=   0 ; % corporate bank risk shock in ESF===> why is this so large?R  %9
 %set_sigma_epsiSH= set_sigma_epsiSF ; % mortgage bank risk shock is ESH %10
 
   set_sigma_epsiWe=data.parameters.sigma_epsiWe;%  [0.00423589602188688]    ; %  entrepreneur net worth shock in EWe %13
  set_sigma_epsiWb=data.parameters.sigma_epsiWb;% banker net worth shock in EWb %11
  set_sigma_epsiRW=data.parameters.sigma_epsiRW;%*** REDUNDANT, RW shock in ERW? what is this?  %12, this seems to be shut off in any case
set_sigma_epsiEC=data.parameters.sigma_epsiEC;
set_sigma_epsiECAB=data.parameters.sigma_epsiECAB;
set_sigma_epsiEbH=data.parameters.sigma_epsiEbH;
set_sigma_epsiEbF=data.parameters.sigma_epsiEbF;
set_sigma_epsiEL=data.parameters.sigma_epsiEL;
  
  set_sigma_e1=data.parameters.sigma_e1; %SS 31
set_sigma_m1=data.parameters.sigma_m1; %SS 32 
%   set_sigma_F=[0.0362988394321432]; %SS 33
%  set_sigma_H=[0.0362992039400083];%SS   variance of iid shock to housing/units of capital/portfolio return for banks F/portfolio return for banks  H resp
% set_sigma_F=[0.04]; 
%set_sigma_H=[0.04];
set_sigma_B=data.parameters.sigma_B;



set_pp  = data.parameters.pp; %inactive transaction cost when recovering funds from failed banks %14
set_betta_s=data.parameters.betta_s; %SS patient HH discount factor %15
set_betta_m= data.parameters.betta_m; %SS impatient HH discount factor %16 
set_hab =data.parameters.hab; %SS habit formation %17



set_alphaa=data.parameters.alphaa; % share of capital in output %22
set_rp=data.parameters.rp; %SS*** loan repayment rate of impatient households %23
set_rpe=data.parameters.rpe; %SS***  loan repayment rate of entrepreneurs %24
set_delta_K= data.parameters.delta_K; %SS depreciation rate of housing & capital shocks %25
set_delta_H= data.parameters.delta_H; % SSsame as above %26
set_mu_m=  data.parameters.mu_m;  %SS     what are these two?  %27
set_mu_e= data.parameters.mu_e ; %SS      28
%set_mu_F=  0.3  ;%SS      corporate bank bankruptcy cost--same as below?  %29
%set_mu_H= 0.3; %SS       mortgage bank bankruptcy cost--different from the paper?  %30
set_mu_B=data.parameters.mu_B;

set_varphi_s=data.parameters.varphi_s; %SS  patient HH preference parameter in utility %35
 set_varphi_m= data.parameters.varphi_m; %SS impatient HH preference parameter in utility %36
set_v_s=  data.parameters.v_s; %SS  patient HH preference parameter in utility  %37==> preference on weight on housing
set_v_m=data.parameters.v_m ; %SS  impatient HH preference parameter in utility %38
set_chi_b= data.parameters.chi_b;%0.15  %SS    CC banker preference parameter in utility %39==>fraction of bankers wealth distributed to savers and/or banker consumption share of wealth
set_chi_e= data.parameters.chi_e; %SS   CC entrepreneur preference parameter in utility  %40
  set_eta =data.parameters.eta; %SS patient HH  inverse frisch elasticity of labor supply  %41
 set_a_e =data.parameters.a_e; % ???  %42
set_a_s=data.parameters.a_s; %SS ??? appears in the budget constraint of savers %43
 set_a_b =data.parameters.a_b; %SS ??? appears in the budget constraint of borrowers %44
set_psi_i= data.parameters.psi_i; % parameter is capital adjustment cost function  %45
set_psi_h=data.parameters.psi_h; % same as above. Why are there two of these?  %46
set_rhoA =data.parameters.rhoA;%47 
set_rhoJ =data.parameters.rhoJ; %48
set_rhoH=  data.parameters.rhoH;
set_rhoK = data.parameters.rhoK; %49
 set_rhoSe =data.parameters.rhoSe; %50
set_rhoSm =data.parameters.rhoSm;  %51
set_rhoSB =  data.parameters.rhoSB;  %53
 set_rhoWb = data.parameters.rhoWb; %54
set_rhoWe =data.parameters.rhoWe;% [0.139485190273843] ;  %55

set_rhoHd =data.parameters.rhoHd; %56
set_rhoHk =  data.parameters.rhoHk; %57
set_rhoRW =data.parameters.rhoRW; %*** shock persistence parameters, same notation as standard deviation %58
set_rhoEC=data.parameters.rhoEC;
set_rhoECAB=data.parameters.rhoECAB;
set_rhoEbH=data.parameters.rhoEbH;
set_rhoEbF=data.parameters.rhoEbF;
set_rhoEL=data.parameters.rhoEL;

set_tau_H=data.parameters.tau_H; %SS*** elasticity of substitution for banks, appears in A1, B1, C1, D1, R_mi R_Fi %63
set_tau_F=data.parameters.tau_F; %SS*** same as above %64
set_phiinf1=data.parameters.phiinf1;  %*** reaction to inflation--why is this so low?  %65
 set_kappa =data.parameters.kappa; %*** interest rate smoothing %66
set_nu=data.parameters.nu; %SS*** penalty cost parameter,  appears in A1, C1 %67
 set_psib =data.parameters.psib; %SS*** penalty cost shape parameter, appears in A1, C1 %68

 %interest stickiness
 set_zeta_m=data.parameters.zeta_m ; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
set_zeta_F=  data.parameters.zeta_F; %SS*** Interest rate stickiness,  appears in 51 FOC for bank business lending, 51 is not in use  %60
 %prudential parameters: 
 %capital requirements
set_Cyphi_H =data.parameters.Cyphi_H; %parameter on capital requirement on mortgage banks %18 ==>pro-cyclical component.
set_Cyphi_F=data.parameters.Cyphi_F; % parameter on capital requirement on corporate banks %20
set_phi_Hs=data.parameters.phi_Hs; % capital adequecy ratio? same equation as above, equation shut off?  %19 --> non-pro cyclical component
set_phis=data.parameters.phis; % is this the capital adequecy ratio?  this is equated to phi, phi is not used anywhere??  is it the same as phi_F and phi_H maybe?  %69
set_phi_Fs =data.parameters.phi_Fs; % capital adequecy ratio? REDUNDANT? appears in the same equation as above, the equation is shut off.  %21
%LTV rule
set_epsilonH1s=data.parameters.epsilonH1s; %***LTV LIMIT?  %61
set_epsilonF1s=data.parameters.epsilonF1s; %*** LTV LIMIT?  =epsilonH  & =epsilonF respectively. 
set_LTVHrule=data.parameters.LTVHrule;%***
set_LTVFrule=data.parameters.LTVFrule;%***
set_gamma_y=  data.parameters.gamma_y;%steady_state(introduced as a free parameter) level of output growth
set_gamma_w=data.parameters.gamma_w;
set_gamma_inve=data.parameters.gamma_inve;
set_gamma_c=data.parameters.gamma_c;
set_gamma_dbe=data.parameters.gamma_dbe;%5.94;
set_gamma_dbm=data.parameters.gamma_dbm;
set_gamma_dq_H=data.parameters.gamma_dq_H;
set_gamma_bspH=data.parameters.gamma_bspH;


set_sigma_epsilon_phi=data.parameters.sigma_epsilon_phi;
set_rho_epsilon_phi=data.parameters.rho_epsilon_phi;
set_def_rate_ss=data.parameters.def_rate_ss;
set_rhoSB=data.parameters.rhoSB;
set_sigma_epsiSB=data.parameters.sigma_epsiSB;
set_omikronH=data.parameters.omikronH;
set_omikronF=data.parameters.omikronF;
% NumberOfParameters = M_.param_nbr;
% for ii = 1:NumberOfParameters
%   paramname = deblank(M_.param_names(ii,:));
%   eval([ paramname ' = M_.params(' int2str(ii) ');']);
% end
set_sigma_epsimarkup_m=data.parameters.sigma_epsimarkup_m;
set_sigma_epsimarkup_F=data.parameters.sigma_epsimarkup_F;
set_rho_markup_m=data.parameters.rho_markup_m;
set_rho_markup_F=data.parameters.rho_markup_F;


all_variables=whos;
tosave = cellfun(@isempty, regexp({all_variables.class}, '^matlab\.(ui|graphics)\.'));


delete LTV1_parameter_values.mat;
save LTV1_parameter_values.mat  -regexp ^(?!(variableToExclude1|variableToExclude2)$).


end