function model_GUI2
diary cmd_output.txt;
set ( 0, 'DefaultFigureColor', [1 1 1] )
data = struct('val',0,'diffMax',1);
data.string=[];
data.endo_name=[];
data.exo_name=[];
data.startDate=datenum('01-01-1998');
data.endDate = datenum('01-12-2016');
data.T=72;%sample length
data.Date=linspace(data.startDate,data.endDate,data.T);

data.parameters.sigma_epsiHd=  0.0093 ;  %  shock to housing depreciation cost in EdH and EdK %1
data.parameters.sigma_epsiHk= 0.0032; % shock to capital depreciation 
data.parameters.sigma_epsiA= 0.0023 ; % productivity shock in A %3
data.parameters.sigma_epsiJ=0.0698 ; % housing preference shock in EJ %4
 data.parameters.sigma_epsiK =   0.0 ; % capital investment shock in EK %5
data.parameters.sigma_epsiH = 0.0 ; % housing investment shock in EH %6
data.parameters.sigma_epsiSe=0.1134; % entrepreneur risk shock in ESe %7
data.parameters.sigma_epsiSm=0.0500; % housing risk shock in ESm %8fd
 data.parameters.sigma_epsiWe= 0.0128 ; %  entrepreneur net worth shock in EWe %13
data.parameters.sigma_epsiWb=0.0074  ;% banker net worth shock in EWb %11
data.parameters.sigma_epsiRW=0;%*** REDUNDANT, RW shock in ERW? what is this?  %12, this seems to be shut off in any case
data.parameters.sigma_epsiEC=0.0137 ;
data.parameters.sigma_epsiECAB=0.0299;
data.parameters.sigma_epsiEbH=0;
data.parameters.sigma_epsiEbF=0;
data.parameters.sigma_epsiEL=0.1;
data.parameters.sigma_epsiSB= 0.0267 ;
  
data.parameters.sigma_e1= 0.1084; %SS 31
data.parameters.sigma_m1=0.0656; %SS 32 
data.parameters.sigma_B=0.080000000000000;


data.parameters.betta_s=0.995; %SS patient HH discount factor %15
data.parameters.betta_m=0.9719; %SS impatient HH discount factor %16 
data.parameters.hab =0.6626; %SS habit formation %17
data.parameters.alphaa=0.3;
data.parameters.pp  = 0.1;% share of capital in output %22
data.parameters.rp=0.01;%SS*** loan repayment rate of impatient households %23
data.parameters.rpe=0.05 ;  %SS***  loan repayment rate of entrepreneurs %24
data.parameters.delta_K= 0.03;  %SS depreciation rate of housing & capital shocks %25
data.parameters.delta_H= 0.01  ; % SSsame as above %26
data.parameters.mu_m=   0.3;  %SS     what are these two?  %27
data.parameters.mu_e= 0.3  ; %SS      28
data.parameters.mu_B=0.3;%set_mu_H=0.3152;
data.parameters.varphi_s=1;%SS  patient HH preference parameter in utility %35
data.parameters.varphi_m= 1;%SS impatient HH preference parameter in utility %36
data.parameters.v_s=  0.25  ;%SS  patient HH preference parameter in utility  %37==> preference on weight on housing
data.parameters.v_m=0.5 ;%SS  impatient HH preference parameter in utility %38
data.parameters.chi_b= 0.15;%SS    CC banker preference parameter in utility %39==>fraction of bankers wealth distributed to savers and/or banker consumption share of wealth
data.parameters.chi_e= 0.1; %SS   CC entrepreneur preference parameter in utility  %40
data.parameters.eta =1;%SS patient HH  inverse frisch elasticity of labor supply  %41
data.parameters.a_e =0; % ???  %42
data.parameters.a_s=0.5; %SS ??? appears in the budget constraint of savers %43
data.parameters.a_b =0.5; %SS ??? appears in the budget constraint of borrowers %44
data.parameters.psi_i=  10.288200000000000 ;% parameter is capital adjustment cost function  %45
data.parameters.psi_h= 12.353199999999999; % same as above. Why are there two of these?  %46

data.parameters.rhoA =0.9801 ;%47 
data.parameters.rhoJ =0.9868; %48
data.parameters.rhoH=   0.0 ;
data.parameters.rhoK = 0.0; %49
data.parameters.rhoSe =0.4171; %50
data.parameters.rhoSm =0.4333;  %51
data.parameters.rhoSB = 0.9738 ;  %52
data.parameters.rhoWb =0.9733 ; %54
data.parameters.rhoWe =0.9188 ;  %55
data.parameters.rhoHd =0.9380 ;  %56
data.parameters.rhoHk =0.8423 ;%57
data.parameters.rhoRW =0; %*** shock persistence parameters, same notation as standard deviation %58
data.parameters.rhoEC=0.8506 ;
data.parameters.rhoECAB=0.2898;
data.parameters.rhoEbH=0;
data.parameters.rhoEbF=0;
data.parameters.rhoEL=0;

data.parameters.tau_H=40; %SS*** elasticity of substitution for banks, appears in A1, B1, C1, D1, R_mi R_Fi %63
data.parameters.tau_F=35; %SS*** same as above %64
data.parameters.phiinf1=1.5;  %*** reaction to inflation--why is this so low?  %65
data.parameters.kappa =0 ; %*** interest rate smoothing %66
data.parameters.nu=0.5; %SS*** penalty cost parameter,  appears in A1, C1 %67
data.parameters.psib =5; %SS*** penalty cost shape parameter, appears in A1, C1 %68

 %interest stickiness
data.parameters.zeta_m=0.691226178352873 ; %SS***  Interest rate stickiness, same as below, appears in  R_m R_F %59
data.parameters.zeta_F=0.447236412060053;  %SS*** Interest rate stickiness,  appears in 51 FOC for bank business lending, 51 is not in use  %60
 %prudential parameters: 
 %capital requirements
data.parameters.Cyphi_H =0; %parameter on capital requirement on mortgage banks %18 ==>pro-cyclical component.
data.parameters.Cyphi_F=0; % parameter on capital requirement on corporate banks %20
data.parameters.phi_Hs=0.11; % capital adequecy ratio? same equation as above, equation shut off?  %19 --> non-pro cyclical component
data.parameters.phis=0.11; % is this the capital adequecy ratio?  this is equated to phi, phi is not used anywhere??  is it the same as phi_F and phi_H maybe?  %69
data.parameters.phi_Fs =0.11; % capital adequecy ratio? REDUNDANT? appears in the same equation as above, the equation is shut off.  %21
%LTV rule
data.parameters.epsilonH1s=0.86; %***LTV LIMIT?  %61
data.parameters.epsilonF1s=0.86; %*** LTV LIMIT?  =epsilonH  & =epsilonF respectively. 
data.parameters.LTVHrule=0;%***
data.parameters.LTVFrule=0;%***
data.parameters.gamma_y=0.206154075688884;%steady_state(introduced as a free parameter) level of output growth
data.parameters.gamma_w=0.770063240360480;
data.parameters.gamma_inve=0.130718781934249;
data.parameters.gamma_c=0.228956168330976;
data.parameters.gamma_dbe=5.926388888888890;%5.94;
data.parameters.gamma_dbm=6.538472222222223;
data.parameters.gamma_dq_H=1.522777777777777;
data.parameters.gamma_bspH=0;

data.parameters.sigma_epsilon_phi=0;
data.parameters.rho_epsilon_phi=0;
data.parameters.def_rate_ss=0;
data.parameters.default_ss=0.1;
data.parameters.omikronH=0.01;
data.parameters.omikronF=0.01;
data.parameters.sigma_epsimarkup_m=  0.0004 ;
data.parameters.sigma_epsimarkup_F= 0.0003;
data.parameters.rho_markup_m=  0.996100000000000 ;
data.parameters.rho_markup_F=  0.993400000000000;
%==================================================
data.parameters.Cyphi_H_counterfactual =0; %parameter on capital requirement on mortgage banks %18 ==>pro-cyclical component.
data.parameters.Cyphi_F_counterfactual =0; % parameter on capital requirement on corporate banks %20
data.parameters.phi_Hs_counterfactual =0.11; % capital adequecy ratio? same equation as above, equation shut off?  %19 --> non-pro cyclical component
data.parameters.phis_counterfactual =0.11; % is this the capital adequecy ratio?  this is equated to phi, phi is not used anywhere??  is it the same as phi_F and phi_H maybe?  %69
data.parameters.phi_Fs_counterfactual  =0.11; % capital adequecy ratio? REDUNDANT? appears in the same equation as above, the equation is shut off.  %21
%LTV rule
data.parameters.epsilonH1s_counterfactual =0.86; %***LTV LIMIT?  %61
data.parameters.epsilonF1s_counterfactual =0.86; %*** LTV LIMIT?  =epsilonH  & =epsilonF respectively. 
data.parameters.LTVHrule_counterfactual =0;%***
data.parameters.LTVFrule_counterfactual =0;%***

%======================================================



f = figure;
    f.Units='normalized';
    f.Position=[0 -0.5 0.8 0.8];
drawnow;

f.Name = 'Sectoral DSGE Model';
tabgp = uitabgroup(f);
tab1 = uitab(tabgp,'Title','Model Setup');
tab2 = uitab(tabgp,'Title','Impulse Responses');
tab3 = uitab(tabgp,'Title','Smoothed Variables');
tab4 = uitab(tabgp,'Title','Historical Variance Decompositions');
 tab5=uitab(tabgp,'Title','Policy Counterfactuals');
 tab6=uitab(tabgp,'Title','Variable and Parameter Definitions');

     ha1 = axes('parent',tab3,'Units','Pixels','Position', [637 162 678 608]);    
     title('Smoothed Variables');
    ha2 = axes('parent',tab2,'Units','Pixels','Position', [637 162 678 608]);    
     title('Impulse Responses');
     ha4 = axes('parent',tab4,'Units','Pixels','Position', [637 162 678 608]);    
     title('Historical Shock Decomposition');
    ha5 = axes('parent',tab5,'Units','Pixels','Position', [637 162 678 250]);    
     title('Counterfactuals: Smoothed Variable');
     
         ha6 = axes('parent',tab5,'Units','Pixels','Position', [637 462 678 250]);    
     title('Counterfactuals: Impulse Response');
     
     uicontrol('parent',tab3,'Style','pushbutton','String','CLEAR FIGURE',...
'Position',  [111 694 189 85],...
'Callback',@clear_smoothedvar_plot);
     uicontrol('parent',tab2,'Style','pushbutton','String','CLEAR FIGURE',...
'Position',  [111 694 189 85],...
'Callback',@clear_irf_plot);
     uicontrol('parent',tab4,'Style','pushbutton','String','CLEAR FIGURE',...
'Position',  [111 694 189 85],...
'Callback',@clear_decomp_plot);
     
 cmd_output= uicontrol('Parent',tab1,'Style','edit',...
    'String','command window output',...
   'Position',[990 100 300 300],'BackgroundColor',[0.5,0.5,0.5]);
cmd_output.String=fileread('cmd_output.txt');
cmd_output.String=cmd_output.String(end-100:end);
     

  uicontrol('parent',tab4,'Style','pushbutton','String','Shock decomposition',...
'Position',  [311 694 189 85],...
'Callback',@compute_shock_decomposition);    
     

uicontrol('parent',tab1,'Style','pushbutton','String','RUN MODEL',...
'Position',  [211 694 189 85],...
'Callback',@run_model_callback);

uicontrol('parent',tab5,'Style','pushbutton','String','RUN COUNTERFACTUAL',...
'Position',  [011 694 189 85],...
'Callback',@run_model_counterfactual_callback);

uicontrol('parent',tab1,'Style','pushbutton','String','RETRIEVE LAST MODEL RUN ',...
'Position',  [511 694 189 85],...
'Callback',@retrieve_last_run);

uicontrol('parent',tab1,'Style','pushbutton','String','SAVE CURRENT FIGURES',...
'Position',  [711 694 189 85],...
'Callback',@save_all_figures);


uicontrol('parent',tab2,'Style','pushbutton','String','Display IRF',...
'Position',[311 694 189 85],...
'Callback',@display_irf);     

uicontrol('parent',tab5,'Style','pushbutton','String','Display IRF',...
'Position',[211 694 89 85],...
'Callback',@display_irf_counterfactual);     

uicontrol('parent',tab5,'Style','pushbutton','String','Display Variable',...
'Position',[311 694 89 85],...
'Callback',@display_simulated_variable_counterfactual);     

uicontrol('parent',tab5,'Style','pushbutton','String','Clear Figures',...
'Position',[411 694 89 85],...
'Callback',@clear_figures_counterfactual);     



uicontrol('parent',tab3,'Style','pushbutton','String','Display Smoothed Variable',...
'Position',[311 694 189 85],...
'Callback',@display_simulated_variable);     


%=======parameter inputs
uicontrol('parent',tab1,'Style','text','String','Sigma e_1','Position',[450 600 40 40]); 
 sigma_e1=uicontrol('parent',tab1,'Style','edit','Callback',@input_sigma_e1,'Position',[490 600 40 40],...
     'String',data.parameters.sigma_e1); 
 
 uicontrol('parent',tab1,'Style','text','String','Sigma m_1','Position',[450 560 40 40]); 
 sigma_m1=uicontrol('parent',tab1,'Style','edit','Callback',@input_sigma_m1,'Position',[490 560 40 40],...
     'String',data.parameters.sigma_m1); 
 
uicontrol('parent',tab1,'Style','text','String','Sigma B','Position',[450 520 40 40]); 
 sigma_B=uicontrol('parent',tab1,'Style','edit','Callback',@input_sigma_B,'Position',[490 520 40 40],...
     'String',data.parameters.sigma_B); 
 

 
 uicontrol('parent',tab1,'Style','text','String','pp','Position',[450 440 40 40]); 
 pp=uicontrol('parent',tab1,'Style','edit','Callback',@input_pp,'Position',[490 440 40 40],...
     'String',data.parameters.pp); 
      
  uicontrol('parent',tab1,'Style','text','String','betta_m','Position',[450 400 40 40]); 
 betta_m=uicontrol('parent',tab1,'Style','edit','Callback',@input_betta_m,'Position',[490 400 40 40],...
     'String',data.parameters.betta_m); 
      
  uicontrol('parent',tab1,'Style','text','String','rp','Position',[450 400 40 40]); 
 rp=uicontrol('parent',tab1,'Style','edit','Callback',@input_rp,'Position',[490 400 40 40],...
     'String',data.parameters.rp); 
     
  uicontrol('parent',tab1,'Style','text','String','rpe','Position',[450 360 40 40]); 
 rpe=uicontrol('parent',tab1,'Style','edit','Callback',@input_rpe,'Position',[490 360 40 40],...
     'String',data.parameters.rpe); 
 
 uicontrol('parent',tab1,'Style','text','String','Cyphi_H','Position',[550 600 40 40]); 
Cyphi_H=uicontrol('parent',tab1,'Style','edit','Callback',@input_Cyphi_H,'Position',[590 600 40 40],...
     'String',data.parameters.Cyphi_H); 
 
 uicontrol('parent',tab1,'Style','text','String','Cyphi_F','Position',[550 560 40 40]); 
Cyphi_F=uicontrol('parent',tab1,'Style','edit','Callback',@input_Cyphi_F,'Position',[590 560  40 40],...
     'String',data.parameters.Cyphi_F); 

 uicontrol('parent',tab1,'Style','text','String','phi_Hs','Position',[550 520  40 40]); 
phi_Hs=uicontrol('parent',tab1,'Style','edit','Callback',@input_phi_Hs,'Position',[590 520  40 40],...
     'String',data.parameters.phi_Hs); 

 uicontrol('parent',tab1,'Style','text','String','phis','Position',[550 480  40 40]); 
phis=uicontrol('parent',tab1,'Style','edit','Callback',@input_phis,'Position',[590 480 40 40],...
     'String',data.parameters.phis); 

 uicontrol('parent',tab1,'Style','text','String','phi_Fs','Position',[550 440 40 40]); 
phi_Fs=uicontrol('parent',tab1,'Style','edit','Callback',@input_phi_Fs,'Position',[590 440 40 40],...
     'String',data.parameters.phi_Fs); 

 uicontrol('parent',tab1,'Style','text','String','epsilonH1s','Position',[550 400 40 40]); 
epsilonH1s=uicontrol('parent',tab1,'Style','edit','Callback',@input_epsilonH1s,'Position',[590 400 40 40],...
     'String',data.parameters.epsilonH1s); 

 uicontrol('parent',tab1,'Style','text','String','epsilonF1s','Position',[550 400 40 40]); 
epsilonF1s=uicontrol('parent',tab1,'Style','edit','Callback',@input_epsilonF1s,'Position',[590 400 40 40],...
     'String',data.parameters.epsilonF1s); 

 uicontrol('parent',tab1,'Style','text','String','LTVHrule','Position',[550 360 40 40]); 
LTVHrule=uicontrol('parent',tab1,'Style','edit','Callback',@input_LTVHrule,'Position',[590 360 40 40],...
     'String',data.parameters.LTVHrule); 

 uicontrol('parent',tab1,'Style','text','String','LTVFrule','Position',[550 320 40 40]); 
LTVFrule=uicontrol('parent',tab1,'Style','edit','Callback',@input_LTVFrule,'Position',[590 320 40 40],...
     'String',data.parameters.LTVFrule); 
%==============================================================================
 uicontrol('parent',tab5,'Style','text','String','Cyphi_H','Position',[500 200 40 40]); 
Cyphi_H_counterfactual=uicontrol('parent',tab5,'Style','edit','Callback',@input_Cyphi_H_counterfactual,'Position',[540 200 40 40],...
     'String',data.parameters.Cyphi_H_counterfactual); 
 
 uicontrol('parent',tab5,'Style','text','String','Cyphi_F','Position',[500 240 40 40]); 
Cyphi_F_counterfactual=uicontrol('parent',tab5,'Style','edit','Callback',@input_Cyphi_F_counterfactual,'Position',[540 240  40 40],...
     'String',data.parameters.Cyphi_F_counterfactual); 

 uicontrol('parent',tab5,'Style','text','String','phi_Hs','Position',[500 280 40 40]); 
phi_Hs_counterfactual=uicontrol('parent',tab5,'Style','edit','Callback',@input_phi_Hs_counterfactual,'Position',[540 280  40 40],...
     'String',data.parameters.phi_Hs_counterfactual); 

 uicontrol('parent',tab5,'Style','text','String','phis','Position',[500 320 40 40]); 
phis_counterfactual=uicontrol('parent',tab5,'Style','edit','Callback',@input_phis_counterfactual,'Position',[540 320 40 40],...
     'String',data.parameters.phis_counterfactual); 

 uicontrol('parent',tab5,'Style','text','String','phi_Fs','Position',[500 360 40 40]); 
phi_Fs_counterfactual=uicontrol('parent',tab5,'Style','edit','Callback',@input_phi_Fs_counterfactual,'Position',[540 360 40 40],...
     'String',data.parameters.phi_Fs_counterfactual); 

 uicontrol('parent',tab5,'Style','text','String','epsilonH1s','Position',[500 400 40 40]); 
epsilonH1s_counterfactual=uicontrol('parent',tab5,'Style','edit','Callback',@input_epsilonH1s_counterfactual,'Position',[540 400 40 40],...
     'String',data.parameters.epsilonH1s_counterfactual); 

 uicontrol('parent',tab5,'Style','text','String','epsilonF1s','Position',[500 440 40 40]); 
epsilonF1s_counterfactual=uicontrol('parent',tab5,'Style','edit','Callback',@input_epsilonF1s_counterfactual,'Position',[540 440 40 40],...
     'String',data.parameters.epsilonF1s_counterfactual); 

 uicontrol('parent',tab5,'Style','text','String','LTVHrule','Position',[500 480 40 40]); 
LTVHrule_counterfactual=uicontrol('parent',tab5,'Style','edit','Callback',@input_LTVHrule_counterfactual,'Position',[540 480 40 40],...
     'String',data.parameters.LTVHrule_counterfactual); 

 uicontrol('parent',tab5,'Style','text','String','LTVFrule','Position',[500 520 40 40]); 
LTVFrule_counterfactual=uicontrol('parent',tab5,'Style','edit','Callback',@input_LTVFrule_counterfactual,'Position',[540 520 40 40],...
     'String',data.parameters.LTVFrule); 
 
 
 
 
  var_names={'b_e','b_m','C','D','def_rate_e','def_rate_m','def_rate_B','R_D','R_F','R_m','R_tilde_F','R_tilde_H','bsp_F' ,'bsp_H',...           
   'phib','dy_data','dq_H_data'};
 %========================================================================
 %===>non-ordered list here might be problematic later...
data.endo_name=uicontrol('style','listbox','min',1,'max',5,...
                 'unit','pix',...
                 'position',[060 000 120 590],...
                 'fontsize',12,...
                 'fontweight','bold',... 
                 'string',{'L'...
                 'L_s'...
                 'L_m'...
                 'H_s'...
                 'H_m'...
                 'Vs'...
                 'Vm'...
                 'Ve'...
                 'Vb'...
                 'b_e',...
                 'b_m'...
                 'C'...
                 'C_s'...
                 'C_m'...
                 'D'...
                 'R_D'...
                 'R_F'...
                 'R_m'...
                 'R_tilde_F'...
                 'R_tilde_H'...
                 'bsp_H'...
                 'bsp_F'...
                 'phib'...
                 'def_rate_m'...
                 'def_rate_e'...
                 'def_rate_B'...
                 'dy_data' ...
                 'dw_data'...
                    'dc_data'...
                       'dinve_data'...
                      'int_rate_HH_2yearLTV75'...
                      'int_rate_business_data'...
                      'bank_rate_data'...
                      'dq_H_data'...
                      'dbm_data'},...
                 'value',1);
uicontrol('Style','text','String','SELECT VARIABLE(S) (max 5)','Position',[40 594 189 40]); 


data.exo_name=uicontrol('style','listbox','min',1,'max',5,...
                 'unit','pix',...
                 'position',[260 000 120 590],...
                 'fontsize',12,...
                 'fontweight','bold',... 
                 'string',{'A' ...
'J' ...
'K' ...
'Se' ...
'Sm' ...
'SB'...
'Wb' ...
'We' ...
'H' ...
'Hd' ...
'Hk' ...
'markup_m'...
'markup_F'...
'EC'...
'ECAB'...
'EL'...
'EbH'...
'EbF'...
'epsilon_phi' ...
},...
                 'value',1);
uicontrol('Style','text','String','SELECT SHOCK   (for impulse responses only)','Position',[240 594 189 40]); 



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

uicontrol('parent',tab5,'Style','text','String','SELECT POLICY PARAMETER  (for counterfactuals)','Position',[440 594 189 40]); 




%=============================================%=============================================
      %=============================================CALLBACK FUNCTIONS
%=============================================%=============================================      
function run_model_callback(~,~)
    f.Pointer='watch';
    drawnow;
data.allvars=run_model(data);
f.Pointer='arrow';
drawnow;
end
%=============================================%=============================================      
function run_model_counterfactual_callback(~,~)
        f.Pointer='watch';
    drawnow;
data.counterfactual=run_model_counterfactual(data);
data.counterfactual_simul=main_counterfactuals(data);
f.Pointer='arrow';
drawnow;
end
%=============================================%=============================================      
function retrieve_last_run(~,~)
        f.Pointer='watch';
    drawnow;
data.allvars=retrieve_model;
f.Pointer='arrow';
drawnow;
end
%=============================================%=============================================      
function save_all_figures(~,~)
        f.Pointer='watch';
    drawnow;
%===================================
   fig1=ha1;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig1,'smoothed_variable','-dpdf');
    close(fig1);
    
f.Pointer='arrow';
drawnow;
end




%==================================================================================================
    function display_irf(~,~)
            f.Pointer='watch';
    drawnow;
        %possibly multiple endogenous variables
        endo_string=get(data.endo_name,'String');
        endo_val=get(data.endo_name,'Value');
        for jj=1:length(endo_val)
        endo_name{jj}=endo_string{endo_val(jj)};
        end
        
        %only one exogenous variable
        exo_string=get(data.exo_name,'String');
        exo_val=get(data.exo_name,'Value');
        exo_name=exo_string{exo_val};
        
       
%         endo_name = get(data.endo_name,{'string','value'});
%         exo_name = get(data.exo_name,{'string','value'});
 axes(ha2);
        for jj=1:length(endo_val)
        string=[endo_name{jj},'_epsi',exo_name];

 if jj==1
        plot(data.allvars.(string),'lineWidth',3,'color','black');
 else
      plot(data.allvars.(string),'lineWidth',3);
 end
 
        hold on;
        end
titlestring=['Impulse Response: Shock ' exo_name];
  title(titlestring);
  legend(endo_name);
  
  f.Pointer='arrow';
drawnow;
    end
%=============================================%=============================================      

    function display_irf_counterfactual(~,~)
            f.Pointer='watch';
    drawnow;
        %possibly multiple endogenous variables
        endo_string=get(data.endo_name,'String');
        endo_val=get(data.endo_name,'Value');
        for jj=1:length(endo_val)
        endo_name{jj}=endo_string{endo_val(jj)};
        end
        
        %only one exogenous variable
        exo_string=get(data.exo_name,'String');
        exo_val=get(data.exo_name,'Value');
        exo_name=exo_string{exo_val};
        
       
%         endo_name = get(data.endo_name,{'string','value'});
%         exo_name = get(data.exo_name,{'string','value'});
 axes(ha6);
        for jj=1:length(endo_val)
        string=[endo_name{jj},'_epsi',exo_name];


        plot(data.allvars.(string),'lineWidth',3,'color','black');
        hold on;
        plot(data.counterfactual.(string),'lineWidth',3,'color','red','LineStyle','--');
% hold on;
        end
titlestring=['Impulse Response: Shock ' exo_name];
  title(titlestring);
  legend('baseline','counterfactual');
  
  
  f.Pointer='arrow';
drawnow;
    end

%=============================================%=============================================      
    function display_simulated_variable(~,~)
            f.Pointer='watch';
    drawnow;
        endo_string=get(data.endo_name,'String');
        endo_val=get(data.endo_name,'Value');
        for jj=1:length(endo_val)
        endo_name{jj}=endo_string{endo_val(jj)};
        end
       
%         endo_name = get(data.endo_name,{'string','value'});
%         exo_name = get(data.exo_name,{'string','value'});
          axes(ha1);
    for jj=1:length(endo_val)
        string=[endo_name{jj}];
       
        
          if jj==1
        plot(data.Date,data.allvars.oo_.SmoothedVariables.(string),'lineWidth',3,'color','black');
          else
               plot(data.Date,data.allvars.oo_.SmoothedVariables.(string),'lineWidth',3);
          end
          
        hold on;
          xlim([data.startDate data.endDate])
  datetick('x','yy','keeplimits');
  titlestring=['Smoothed Variables'];
  title(titlestring);
    legend(endo_name);
    end

    
    f.Pointer='arrow';
drawnow;
    end
%=============================================%=============================================      
    function display_simulated_variable_counterfactual(~,~)
            f.Pointer='watch';
    drawnow;
        endo_string=get(data.endo_name,'String');
        endo_val=get(data.endo_name,'Value');
        for jj=1:length(endo_val)
        endo_name{jj}=endo_string{endo_val(jj)};
        end
       

          axes(ha5);
          
   field_names=cellstr(data.allvars.M_.endo_names);
 index=1;
  for jj=1:length(field_names);
  
            if true(strcmp(field_names(jj),endo_name(index)))
var_indices(index)=jj;
                     if index<length(endo_name)
                    index=index+1;
                     end

            end
      
  end        
          
 %this won't work in multiple variables are selected
 

 
 
    plot(data.Date,data.counterfactual_simul.baseline(var_indices(1),:),'lineWidth',3,'color','black');
    hold on;
    plot(data.Date,data.counterfactual_simul.counterfactual(var_indices(1),:),'lineWidth',3,'color','red','LineStyle','--');
          xlim([data.startDate data.endDate])
  datetick('x','yy','keeplimits');
  titlestring=['Smoothed Variables'];
  title(titlestring);
    legend('baseline','counterfactual');
          
%     for jj=1:length(endo_val)
%         string=[endo_name{jj}];
%        
%         
%           if jj==1
%         plot(data.Date,data.counterfactuals.baseline(:,jj),'lineWidth',3,'color','black');
%         hold on;
%         plot(data.Date,data.counterfactuals.counterfactual(:,jj),'lineWidth',3,'color','black');
%           else
%                plot(data.Date,data.allvars.oo_.SmoothedVariables.(string),'lineWidth',3);
%           end
   
f.Pointer='arrow';
drawnow;
    end

%================================================================================================

function compute_shock_decomposition(~,~)
         f.Pointer='watch';
    drawnow;

    endo_string=get(data.endo_name,'String');
        endo_val=get(data.endo_name,'Value');
        for jj=1:length(endo_val)
        endo_name{jj}=endo_string{endo_val(jj)};
        end   
        


M_=data.allvars.M_;
oo_=data.allvars.oo_;
options_=data.allvars.options_;
shock_decomp=oo_.shock_decomposition;

 field_names=cellstr(M_.endo_names);
 index=1;
  for jj=1:length(field_names);
  
            if true(strcmp(field_names(jj),endo_name(index)))
var_indices(index)=jj;
                     if index<length(endo_name)
                    index=index+1;
                     end

            end
      
  end

data_for_plot=squeeze(shock_decomp(var_indices,:,:))
data_for_plot=data_for_plot';
axes(ha4);
bar(data_for_plot);



  f.Pointer='arrow';
drawnow;
end
%=============================================%=============================================      

function plot_counterfactual_irf(~,~)
        f.Pointer='watch';
    drawnow;
         endo_string=get(data.endo_name,'String');
        endo_val=get(data.endo_name,'Value');
        for jj=1:length(endo_val)
        endo_name{jj}=endo_string{endo_val(jj)};
        end   
        
        
          axes(ha5);
    for jj=1:length(endo_val)
        string=[endo_name{jj}];
       
        
          if jj==1
        plot(data.Date,data.allvars.oo_.SmoothedVariables.(string),'lineWidth',3,'color','black');
        hold on;
        plot(data.Date,data.counterfactual.oo_.SmoothedVariables.(string),'lineWidth',3,'color','red','LineStyle','--');
          else
               plot(data.Date,data.allvars.oo_.SmoothedVariables.(string),'lineWidth',3);
               hold on;
               plot(data.Date,data.counterfactual.oo_.SmoothedVariables.(string),'lineWidth',3,'LineStyle','--');
          end
          
        hold on;
          xlim([data.startDate data.endDate])
  datetick('x','yy','keeplimits');
  titlestring=['Smoothed Variables'];
  title(titlestring);
    legend(endo_name);
    end        
        
        
        
f.Pointer='arrow';
drawnow;
end


%=============================================%=============================================      


    function clear_smoothedvar_plot(~,~)
            f.Pointer='watch';
    drawnow;
        cla(ha1, 'reset');
        ha1.Parent=tab3;
        ha1.Units='Pixels';
        ha1.Position=[537 162 878 608];
     title('Smoothed Variables');
     legend('');
     f.Pointer='arrow';
drawnow;
    end
%=============================================%=============================================      

    function clear_irf_plot(~,~)
            f.Pointer='watch';
    drawnow;
         cla(ha2, 'reset');
        ha2.Parent=tab2;
        ha2.Units='Pixels';
        ha2.Position=[537 162 878 608];
     title('Impulse Responses');
     legend('');
     f.Pointer='arrow';
drawnow;
    end
%=============================================%=============================================      
    function clear_decomp_plot(~,~)
            f.Pointer='watch';
    drawnow;
         cla(ha4, 'reset');
        ha4.Parent=tab4;
        ha4.Units='Pixels';
        ha4.Position=[537 162 878 608];
     title('Compute Shock Decomposition');
     legend('');
     f.Pointer='arrow';
drawnow;
    end
%=============================================%=============================================      

    function clear_figures_counterfactual(~,~)
      f.Pointer='watch';
    drawnow;
         cla(ha5,'reset');
        ha5.Parent=tab5;
        ha5.Units='Pixels';
         ha5.Position=[637 162 678 250];
        title('Counterfactuals: Smoothed Variable');
            legend('');
     
     
        cla(ha6, 'reset');
        ha6.Parent=tab5;
        ha6.Units='Pixels';
        ha6.Position=[637 462 678 250];
         title('Counterfactuals: Impulse Response');
        legend('');
        f.Pointer='arrow';
drawnow;
    end



%=============================================%=============================================      
%=====================
%=========input functions============
%==leave these at the bottom plz
%=============================================%=============================================      


 function input_sigma_e1(~,~)
data.parameters.sigma_e1=str2double(get(sigma_e1,'String'));
 end

 function input_sigma_m1(~,~)
data.parameters.sigma_m1=str2double(get(sigma_m1,'String'));
 end

 function input_sigma_B(~,~)
data.parameters.sigma_B=str2double(get(sigma_B,'String'));
 end



 function input_pp(~,~)
data.parameters.pp=str2double(get(pp,'String'));
 end

function input_betta_m(~,~)
data.parameters.betta_m=str2double(get(betta_m,'String'));
 end

 function input_rp(~,~)
data.parameters.rp=str2double(get(rp,'String'));
 end

 function input_rpe(~,~)
data.parameters.rpe=str2double(get(rpe,'String'));
 end

 function input_Cyphi_H(~,~)
data.parameters.Cyphi_H=str2double(get(Cyphi_H,'String'));
 end

function input_Cyphi_F(~,~)
data.parameters.Cyphi_F=str2double(get(Cyphi_F,'String'));
end

function input_phi_Hs(~,~)
data.parameters.phi_Hs=str2double(get(phi_Hs,'String'));
end

function input_phis(~,~)
data.parameters.phis=str2double(get(phis,'String'));
end

function input_phi_Fs(~,~)
data.parameters.phi_Fs=str2double(get(phi_Fs,'String'));
end

function input_epsilonH1s(~,~)
data.parameters.epsilonH1s=str2double(get(epsilonH1s,'String'));
end

function input_epsilonF1s(~,~)
data.parameters.epsilonF1s=str2double(get(epsilonF1s,'String'));
end

function input_LTVHrule(~,~)
data.parameters.LTVHrule=str2double(get(LTVHrule,'String'));
end

function input_LTVFrule(~,~)
data.parameters.LTVFrule=str2double(get(LTVFrule,'String'));
end

%==================================================
%---------------input functions for counterfactual parameters--------------------
function input_Cyphi_H_counterfactual(~,~)
data.parameters.Cyphi_H_counterfactual=str2double(get(Cyphi_H_counterfactual,'String'));
 end

function input_Cyphi_F_counterfactual(~,~)
data.parameters.Cyphi_F_counterfactual=str2double(get(Cyphi_F_counterfactual,'String'));
end

function input_phi_Hs_counterfactual(~,~)
data.parameters.phi_Hs_counterfactual=str2double(get(phi_Hs_counterfactual,'String'));
end

function input_phis_counterfactual(~,~)
data.parameters.phis_counterfactual=str2double(get(phis_counterfactual,'String'));
end

function input_phi_Fs_counterfactual(~,~)
data.parameters.phi_Fs_counterfactual=str2double(get(phi_Fs_counterfactual,'String'));
end

function input_epsilonH1s_counterfactual(~,~)
data.parameters.epsilonH1s_counterfactual=str2double(get(epsilonH1s_counterfactual,'String'));
end

function input_epsilonF1s_counterfactual(~,~)
data.parameters.epsilonF1s_counterfactual=str2double(get(epsilonF1s_counterfactual,'String'));
end

function input_LTVHrule_counterfactual(~,~)
data.parameters.LTVHrule_counterfactual=str2double(get(LTVHrule_counterfactual,'String'));
end

function input_LTVFrule_counterfactual(~,~)
data.parameters.LTVFrule_counterfactual=str2double(get(LTVFrule_counterfactual,'String'));
end

end