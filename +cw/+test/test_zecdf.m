function test_zecdf
    d = exprnd(30,3000,1);
    
    [CDF_z,t_z] = zecdf(d);
    [CDF_m,t_m] = ecdf(d);
    
    figure(19783)
    clf
    
        subplot(2,1,1)
            hold all
            
            hist(d,50)
            
        subplot(2,1,2)
            hold all
            
            plot(t_m,CDF_m,'r-','LineWidth',3)
            plot(t_z,CDF_z,'g--','LineWidth',3)
            
            set(gca, 'XScale', 'log')
end

