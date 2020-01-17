%M=mortgage bank, E=entrepreneur, H=household, F=corporate bank
%***starred variables are not present in the original 3-D model
%CC variables denote calibrates ones in CC
set_sigma_epsiHd=   2.0856  ;  %  shock to housing depreciation cost in EdH and EdK %1
set_sigma_epsiHk= 0.7480 ; % shock to capital depreciation 

 set_sigma_epsiA=    1.6160  ; % productivity shock in A %3
 set_sigma_epsiJ=  0.0796   ; % housing preference shock in EJ %4
 
 set_sigma_epsiK =  0.0799  ; % capital investment shock in EK %5
set_sigma_epsiH =0.0789  ; % housing investment shock in EH %6

set_sigma_epsiSe=0.0447     ; % entrepreneur risk shock in ESe %7
set_sigma_epsiSm= 0.0449    ; % housing risk shock in ESm %8fd
  set_sigma_epsiSF=   0.0506     ; % corporate bank risk shock in ESF===> why is this so large?R  %9
 set_sigma_epsiSH= set_sigma_epsiSF ; % mortgage bank risk shock is ESH %10
 
   set_sigma_epsiWe= 0.0448      ; %  entrepreneur net worth shock in EWe %13
  set_sigma_epsiWb=0.0450 ;% banker net worth shock in EWb %11
  set_sigma_epsiRW=0;%*** REDUNDANT, RW shock in ERW? what is this?  %12, this seems to be shut off in any case

  
  set_sigma_e1=0.03; %31
set_sigma_m1=0.03; %32 
set_sigma_F=0.03; %33
set_sigma_H=set_sigma_F;%  variance of iid shock to housing/units of capital/portfolio return for banks F/portfolio return for banks  H resp
  
set_pp  = 0 ; % transaction cost when recovering funds from failed banks %14
set_betta_s=0.995; % patient HH discount factor %15
set_betta_m=0.9598; % impatient HH discount factor %16 
set_hab =  0.7996   ; % habit formation %17



set_alphaa= 0.3; % share of capital in output %22
set_rp=0.03; %*** loan repayment rate of impatient households %23
set_rpe= 0.03; %***  loan repayment rate of entrepreneurs %24
set_delta_K=0.03; % depreciation rate of housing & capital shocks %25
set_delta_H= 0.03; % same as above %26
set_mu_m=  0.2985;  %what are these two?  %27
set_mu_e=0.2986  ; %28
set_mu_F=0.2984  ;%corporate bank bankruptcy cost--same as below?  %29
set_mu_H= set_mu_F; %mortgage bank bankruptcy cost--different from the paper?  %30

set_varphi_s=    1.0037 ; %  patient HH preference parameter in utility %35
 set_varphi_m=0.9974; % impatient HH preference parameter in utility %36
set_v_s=  0.0981 ; %  patient HH preference parameter in utility  %37==> preference on weight on housing
set_v_m=0.0039; %  impatient HH preference parameter in utility %38
set_chi_b= 0.2194; %CC banker preference parameter in utility %39==>fraction of bankers wealth distributed to savers and/or banker consumption share of wealth
set_chi_e= 0.3906; %CC entrepreneur preference parameter in utility  %40
  set_eta =0.2; % patient HH  inverse frisch elasticity of labor supply  %41
 set_a_e =0; % ???  %42
set_a_s=0.5; % ??? appears in the budget constraint of savers %43
 set_a_b =0.5; % ??? appears in the budget constraint of borrowers %44
set_psi_i=  2.0217 ; % parameter is capital adjustment cost function  %45
set_psi_h=2.0519  ; % same as above. Why are there two of these?  %46
set_rhoA =0.9;%47 
set_rhoJ =0.6974  ; %48
set_rhoH=   0.6974 ;
set_rhoK = 0.6975 ; %49
 set_rhoSe =    0.5000    ; %50
set_rhoSm =   0.5000  ;  %51
set_rhoSF =    0.5000  ;  %52
set_rhoSH =  set_rhoSF ;  %53
 set_rhoWb = 0.5000; %54
set_rhoWe =  0.5099 ;  %55

set_rhoHd =  0.9946  ; %56
set_rhoHk =  0.7470 ; %57
set_rhoRW =0; %*** shock persistence parameters, same notation as standard deviation %58


set_tau=50; %*** elasticity of substitution for banks, appears in A1, B1, C1, D1, R_mi R_Fi %63
set_taue=50; %*** same as above %64
set_phiinf1=1.5;  %*** reaction to inflation--why is this so low?  %65
 set_kappa =0 ; %*** interest rate smoothing %66
set_nu=0.3; %*** penalty cost parameter,  appears in A1, C1 %67
 set_psib =10; %*** penalty cost shape parameter, appears in A1, C1 %68

 %interest stickiness
 set_zeta1=0.02; %***  Interest rate stickiness, same as below, appears in  R_m R_F %59
set_zetae=set_zeta1; %*** Interest rate stickiness,  appears in 51 FOC for bank business lending, 51 is not in use  %60
 %prudential parameters: 
 %capital requirements
set_Cyphi_H =0; %parameter on capital requirement on mortgage banks %18 ==>pro-cyclical component.
set_Cyphi_F=set_Cyphi_H; % parameter on capital requirement on corporate banks %20
set_phi_Hs=0.11; % capital adequecy ratio? same equation as above, equation shut off?  %19 --> non-pro cyclical component
set_phis=set_phi_Hs; % is this the capital adequecy ratio?  this is equated to phi, phi is not used anywhere??  is it the same as phi_F and phi_H maybe?  %69
set_phi_Fs =set_phi_Hs; % capital adequecy ratio? REDUNDANT? appears in the same equation as above, the equation is shut off.  %21
%LTV rule
set_epsilonH1s=0.86; %***LTV LIMIT?  %61
set_epsilonF1s=0.86; %*** LTV LIMIT?  =epsilonH  & =epsilonF respectively. 
set_LTVHrule=0;%***
set_LTVFrule=0;%***
set_gamma_y= 0.3292 ;%steady_state(introduced as a free parameter) level of output growth
set_gamma_db=5.94;
delete LTV1_parameter_values.mat;
save LTV1_parameter_values.mat;
