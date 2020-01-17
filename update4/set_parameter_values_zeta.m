function []=set_parameter_values_zeta(set_zeta_m,set_zeta_F)
set_sigma_epsiHd=    0.0093 ;  %  shock to housing depreciation cost in EdH and EdK %1
set_sigma_epsiHk= 0.0032 ; % shock to capital depreciation 
set_sigma_epsiA=   0.0023; % productivity shock in A %3
 set_sigma_epsiJ= 0.5741   ; % housing preference shock in EJ %4
set_sigma_epsiK =   0.0 ; % capital investment shock in EK %5
set_sigma_epsiH = 0.0 ; % housing investment shock in EH %6
set_sigma_epsiSe= 0.0492  ; % entrepreneur risk shock in ESe %7
set_sigma_epsiSm=  0.0481 ; % housing risk shock in ESm %8fd
 set_sigma_epsiWe=0.0069 ;%  [0.00423589602188688]    ; %  entrepreneur net worth shock in EWe %13
  set_sigma_epsiWb=0 ;% banker net worth shock in EWb %11
  set_sigma_epsiRW=0;%*** REDUNDANT, RW shock in ERW? what is this?  %12, this seems to be shut off in any case
set_sigma_epsiEC=0.0155 ;
set_sigma_epsiECAB=0.0698;
set_sigma_epsiEbH=0;
set_sigma_epsiEbF=0;
set_sigma_epsiEL=0;
set_sigma_e1=  0.1084 ; 
set_sigma_m1=  0.0656 ;
set_sigma_B=   0.0679;
set_sigma_epsimarkup_m=  0.0004;
set_sigma_epsimarkup_F= 0.0005;
set_sigma_epsiSB=0.0267;
set_sigma_epsilon_phi=0;


set_rhoA = 0.9775  ;%47 
set_rhoJ =  0.9227; %48
set_rhoH=   0.0 ;
set_rhoK = 0.0; %49
 set_rhoSe = 0.5619 ; %50
set_rhoSm =0.5427  ;  %51
set_rhoSF =    0 ;  %52
set_rhoSH =  set_rhoSF ;  %53
 set_rhoWb =0; %54
set_rhoWe =  0.2284;% [0.139485190273843] ;  %55
set_rhoHd =0.9300 ; %56
set_rhoHk = 0.9919 ; %57
set_rhoRW =0; %*** shock persistence parameters, same notation as standard deviation %58
set_rhoEC= 0.9139;
set_rhoECAB=0.5260  ;
set_rhoEbH=0;
set_rhoEbF=0;
set_rhoEL=0;
set_rho_markup_m=0.9961;
set_rho_markup_F=0.9934;
set_rho_epsilon_phi=0;
set_rhoSB=0.0570;

set_pp  = 0.1; %inactive transaction cost when recovering funds from failed banks %14
set_betta_s=0.995; %SS patient HH discount factor %15
set_betta_m= 0.9719; %SS impatient HH discount factor %16 
set_hab =0.6626; %SS habit formation %17
set_alphaa= 0.3 ; % share of capital in output %22
set_rp=0.0100; %SS*** loan repayment rate of impatient households %23
set_rpe=  0.05; %SS***  loan repayment rate of entrepreneurs %24
set_delta_K= 0.03; %SS depreciation rate of housing & capital shocks %25
set_delta_H= 0.01 ; % SSsame as above %26
set_mu_m=   0.3;  %SS     what are these two?  %27
set_mu_e= 0.3  ; %SS      28
set_mu_B=0.3;
set_varphi_s=1; %SS  patient HH preference parameter in utility %35
 set_varphi_m=  1; %SS impatient HH preference parameter in utility %36
set_v_s=    0.25  ; %SS  patient HH preference parameter in utility  %37==> preference on weight on housing
set_v_m=0.5 ; %SS  impatient HH preference parameter in utility %38
set_chi_b= 0.15;%0.15  %SS    CC banker preference parameter in utility %39==>fraction of bankers wealth distributed to savers and/or banker consumption share of wealth
set_chi_e= 0.1; %SS   CC entrepreneur preference parameter in utility  %40
  set_eta =1; %SS patient HH  inverse frisch elasticity of labor supply  %41
 set_a_e =0; % ???  %42
set_a_s=1/2; %SS ??? appears in the budget constraint of savers %43
 set_a_b =1/2; %SS ??? appears in the budget constraint of borrowers %44
set_psi_i=10.2882 ; % parameter is capital adjustment cost function  %45
set_psi_h=12.3532; % same as above. Why are there two of these?  %46


%STICKINESS & CAPITAL REQ. RELATED STUFF
set_tau_H=40; %SS*** elasticity of substitution for banks, appears in A1, B1, C1, D1, R_mi R_Fi %63
set_tau_F=35; %SS*** same as above %64
set_phiinf1=1.5;  %*** reaction to inflation--why is this so low?  %65
 set_kappa =0 ; %*** interest rate smoothing %66
set_nu=0.5; %SS*** penalty cost parameter,  appears in A1, C1 %67
 set_psib =05; %SS*** penalty cost shape parameter, appears in A1, C1 %68

 %interest stickiness
%set_zeta_m=[0.691226178352873]; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
%set_zeta_F=[0.447236412060053]; %SS*** Interest rate stickiness,  appears in 51 FOC for bank business lending, 51 is not in use  %60
 %prudential parameters: 
 %capital requirements


%HISTORICAL MEAN PARAMETERS FOR MEASUREMENT EQUATIONS
set_gamma_y=0.206154075688884;%steady_state(introduced as a free parameter) level of output growth
set_gamma_w=0.770063240360480;
set_gamma_inve=0.130718781934249;
set_gamma_c=0.228956168330976;
set_gamma_dbe=5.926388888888890;%5.94;
set_gamma_dbm=6.538472222222223;
set_gamma_dq_H= 1.522777777777777;
set_gamma_bspH=0;
set_def_rate_ss=0;





%===============PRUDENTIAL POLICY PARAMETERS
set_Cyphi_H =0; %parameter on capital requirement on mortgage banks %18 ==>pro-cyclical component.
set_Cyphi_F=0; % parameter on capital requirement on corporate banks %20
set_phi_Hs=0.11; % capital adequecy ratio? same equation as above, equation shut off?  %19 --> non-pro cyclical component
set_phi_Fs =0.11; % capital adequecy ratio? REDUNDANT? appears in the same equation as above, the equation is shut off.  %21
set_phis=0.11; % is this the capital adequecy ratio?  this is equated to phi, phi is not used anywhere??  is it the same as phi_F and phi_H maybe?  %69
%LTV rule
set_epsilonH1s=0.86; %***LTV LIMIT?  %61
set_epsilonF1s=0.86; %*** LTV LIMIT?  =epsilonH  & =epsilonF respectively. 
set_LTVHrule=0;%***
set_LTVFrule=0;%***
set_omikronH=0.01;
set_omikronF=0.01;


delete LTV1_parameter_values.mat;
save LTV1_parameter_values.mat;

end