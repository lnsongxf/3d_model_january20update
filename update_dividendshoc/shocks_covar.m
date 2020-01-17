smoothed_shocks=[oo_.SmoothedShocks.epsiA',
    oo_.SmoothedShocks.epsiJ',
     oo_.SmoothedShocks.epsiSe',
      oo_.SmoothedShocks.epsiSB',
       oo_.SmoothedShocks.epsiWe',
        oo_.SmoothedShocks.epsiWb',
         oo_.SmoothedShocks.epsimarkup_m',
          oo_.SmoothedShocks.epsimarkup_F',
           oo_.SmoothedShocks.epsiEC',
            oo_.SmoothedShocks.epsiECAB',
             oo_.SmoothedShocks.epsiEbF']';
         
         
  [shocks_covar_ pval_]=corr(smoothed_shocks)       