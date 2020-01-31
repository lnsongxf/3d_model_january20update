LTV1_steadystate
LTV1_steadystate
dynare LTV1.mod
dynare LTV1.mod
load('LTV1_results.mat', 'oo_')
irf1=oo_.endo_simul;
save irf1
load irf1

LTV1INT_steadystate
LTV1INT_steadystate
dynare LTV1INT.mod
dynare LTV1INT.mod
load('LTV1INT_results.mat', 'oo_')
irf2=oo_.endo_simul;
save irf2
load irf2

load irf1



b_m1=irf1(7,:);
db_m1=irf1(7,:)-irf1(7,1);
pb_m1=100*((irf1(7,:)-irf1(7,1)))/irf1(7,1);
R_m1=irf1(62,:);
dR_m1=irf1(62,:)-irf1(62,1);
b_e1=irf1(1,:);
db_e1=irf1(1,:)-irf1(1,1);
dR_F1=irf1(58,:)-irf1(58,1);
pb_e1=100*((irf1(1,:)-irf1(1,1)))/irf1(1,1);
b_tot1=irf1(162,:);
q_H1=irf1(55,:);
q_K1=irf1(56,:);
H_m1=irf1(39,:);
pH_m1=100*((irf1(39,:)-irf1(39,1)))/irf1(39,1);
K1=irf1(42,:);
pK1=100*((irf1(42,:)-irf1(42,1)))/irf1(42,1);



b_m2=irf2(7,:);
db_m2=irf2(7,:)-irf2(7,1);
pb_m2=100*((irf2(7,:)-irf2(7,1)))/irf2(7,1);
R_m2=irf2(62,:);
dR_m2=irf2(62,:)-irf2(62,1);
b_e2=irf2(1,:);
dR_F2=irf2(58,:)-irf2(58,1);
db_e2=irf2(1,:)-irf2(1,1);
pb_e2=100*((irf2(1,:)-irf2(1,1)))/irf2(1,1);
b_tot2=irf2(162,:);
q_H2=irf2(55,:);
q_K2=irf2(56,:);
H_m2=irf2(39,:);
pH_m2=100*((irf2(39,:)-irf2(39,1)))/irf2(39,1);
K2=irf2(42,:);
pK2=100*((irf2(42,:)-irf2(42,1)))/irf2(42,1);




  figure;
          HOR=1:1:99;
 graph23 =plot(HOR,dR_F1(2:100),HOR,dR_F2(2:100)) ;
        set(graph23,'Linewidth',3);
        title(['Business Interest Rates'] )
        legend(['dR_F1', 'dR_F2'], {'\xi=0.02', '\xi=0.45'})


figure;
HOR=1:1:100;
 graph22 =plot(HOR,b_e1(1:100),HOR,b_e2(1:100)) ;
        set(graph22,'Linewidth',3);
        title(['Business Loans (SS Level)'] )
         legend(['b_e1', 'b_e2'], {'\xi=0.02', '\xi=0.45'})

figure;
HOR=1:1:100;
 graph1 =plot(HOR,b_m1(1:100),HOR,b_m2(1:100)) ;
        set(graph1,'Linewidth',3);
        title(['Mortgage Loans (SS levels)'] )
         legend(['b_m1', 'b_m2'],{'\xi=0.02', '\xi=0.45'})
  figure;
 HOR=1:1:100;
 graph9 =plot(HOR,db_m1(1:100),HOR,db_m2(1:100)) ;
        set(graph9,'Linewidth',3);
        title(['Mortgage Loans'] )
        legend(['db_m1', 'db_m2'],{'\xi=0.02', '\xi=0.45'})  
 
        figure;
          HOR=1:1:99;
 graph10 =plot(HOR,pb_m1(2:100),HOR,pb_m2(2:100)) ;
        set(graph10,'Linewidth',3);
        title(['Mortgage Loans'] )
        legend(['db_m1', 'db_m2'],{'\xi=0.02', '\xi=0.45'}) 
        
        
        
     
        
        
              
         figure;    
   graph13 =plot(b_m1(1:100)) ;
        set(graph13,'Linewidth',3);
        title(['Mortgage Loans'] )
        
     
     figure;
          HOR=1:1:99;
 graph14 =plot(HOR,R_m1(2:100),HOR,R_m2(2:100)) ;
        set(graph14,'Linewidth',3);
        title(['Mortgage Interest Rates'] )
        legend(['R_m1', 'R_m2'],{'\xi=0.02', '\xi=0.45'})   
        
     
           figure;
          HOR=1:1:99;
 graph15 =plot(HOR,dR_m1(2:100),HOR,dR_m2(2:100)) ;
        set(graph15,'Linewidth',3);
        title(['Mortgage Interest Rates'] )
        legend(['dR_m1', 'dR_m2'],{'\xi=0.02', '\xi=0.45'})
        
        
        
        
    figure;    
   graph2 =plot(R_m1(1:100)) ;
        set(graph2,'Linewidth',3);
        title(['Interest Rates'] )
              
        figure;
        graph3 =plot(b_e1(1:100)) ;
        set(graph3,'Linewidth',3);
        title(['Business Loans'] )
        
          figure;
          HOR=1:1:100;
 graph16 =plot(HOR,db_e1(1:100),HOR,db_e2(1:100)) ;
        set(graph16,'Linewidth',3);
        title(['Business Loan'] )
        legend(['db_e1', 'db_e2'], {'\xi=0.02', '\xi=0.45'})
        
          figure;
          HOR=1:1:99;
 graph16 =plot(HOR,pb_e1(2:100),HOR,pb_e2(2:100)) ;
        set(graph16,'Linewidth',3);
        title(['Business Loan'] )
        legend(['pb_e1', 'pb_e2'],{'\xi=0.02', '\xi=0.45'})
        
        
        figure;
        graph4 =plot(b_tot1(1:100)) ;
        set(graph4,'Linewidth',3);
        title(['Total Loans'] )
        
        figure;
          graph5 =plot(q_H1(1:100)) ;
        set(graph5,'Linewidth',3);
        title(['House price'] )
        figure;
          graph6 =plot(q_K1(1:100)) ;
        set(graph6,'Linewidth',3);
        title(['Capital price'] )
        
        figure;
          HOR=1:1:100;
 graph17 =plot(HOR,q_H1(1:100),HOR,q_H2(1:100)) ;
        set(graph17,'Linewidth',3);
        title(['Housing Price'] )
        legend(['q_H1', 'q_H2'],{'\xi=0.02', '\xi=0.45'})
        
         figure;
          HOR=1:1:100;
 graph18 =plot(HOR,q_K1(1:100),HOR,q_K2(1:100)) ;
        set(graph18,'Linewidth',3);
        title(['Capital Price'] )
        legend(['q_K1', 'q_K2'],{'\xi=0.02', '\xi=0.45'})
        
       
         figure;
          HOR=1:1:100;
 graph19 =plot(HOR,H_m1(1:100),HOR,H_m2(1:100)) ;
        set(graph19,'Linewidth',3);
        title(['Housing Stock'] )
        legend(['H_m1', 'H_m2'],{'\xi=0.02', '\xi=0.45'})
        
       
        
        
        
        figure;
          HOR=1:1:99;
 graph20 =plot(HOR,pH_m1(2:100),HOR,pH_m2(2:100)) ;
        set(graph20,'Linewidth',3);
        title(['Housing Stock'] )
        legend(['pH_m1', 'pH_m2'],{'\xi=0.02', '\xi=0.45'})
        
          figure;
          HOR=1:1:99;
 graph21 =plot(HOR,pK1(2:100),HOR,pK2(2:100)) ;
        set(graph21,'Linewidth',3);
        title(['Capital Stock'] )
        legend(['pK1', 'pK2'],{'\xi=0.02', '\xi=0.45'})
        
        figure;
          graph7 =plot(H_m1(1:100)) ;
        set(graph7,'Linewidth',3);
        title(['Housing stock'] )
        figure;
         graph8 =plot(K1(1:100)) ;
        set(graph8,'Linewidth',3);
        title(['Capital'] )
        