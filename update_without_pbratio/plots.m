 figure;subplot(4,1,1);plot(def_rate_e,'color','black');title('def rate e ');
  subplot(4,1,2);plot(def_rate_m,'color','black');title('def rate m ');
   subplot(4,1,3);plot(def_rate_H,'color','black');title('def rate H');
   subplot(4,1,4);plot(def_rate_F,'color','black');title('def rate F');
    
   
    figure;subplot(2,1,1);plot(b_to_Y_gap_data);subplot(2,1,2);autocorr(b_to_Y_gap_data);
figure;subplot(2,1,1);plot(diff(b_to_Y_gap_data));subplot(2,1,2);autocorr(diff(b_to_Y_gap_data));