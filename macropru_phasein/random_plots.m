%run after LTV1.mod
% figure;
% subplot(2,2,1);
% plot(oo_.SmoothedVariables.A1);title('A1');
% subplot(2,2,2);
% plot(oo_.SmoothedVariables.B1);title('B1');
% subplot(2,2,3);
% plot(oo_.SmoothedVariables.C1);title('C1');
% subplot(2,2,4);
% plot(oo_.SmoothedVariables.D1);title('D1');

figure;
subplot(2,2,1);
plot(oo_.SmoothedVariables.R_mi);title('R_{mi}');
subplot(2,2,2);
plot(oo_.SmoothedVariables.R_m);title('R_m');
subplot(2,2,3);
plot(oo_.SmoothedVariables.R_Fi);title('R_{Fi}');
subplot(2,2,4);
plot(oo_.SmoothedVariables.R_F);title('R_F');

figure;
plot(oo_.SmoothedVariables.phib);title('Capital ratio');

figure;
subplot(5,2,1);
plot(oo_.SmoothedVariables.def_rate_e);title('def rate e');
subplot(5,2,2);
plot(oo_.SmoothedVariables.def_rate_m);title('def rate m');
subplot(5,2,3);
plot(oo_.SmoothedVariables.def_rate_B);title('def rate B');
subplot(5,2,4);
plot(oo_.SmoothedVariables.D);title('Deposits');
subplot(5,2,5);
plot(oo_.SmoothedVariables.b_e);title('Business Borrowing');
subplot(5,2,6);
plot(oo_.SmoothedVariables.b_m);title('Household Borrowing');
subplot(5,2,7);
plot(oo_.SmoothedVariables.dbe_data);title('Growth of Business Borrowing');
subplot(5,2,8);
plot(oo_.SmoothedVariables.dbm_data);title('Growth of Household Borrowing');
subplot(5,2,9);
plot(oo_.SmoothedVariables.R_tilde_F);title('Return on Business Portfolio');
subplot(5,2,10);
plot(oo_.SmoothedVariables.R_tilde_H);title('Return on Housing Portfolio')