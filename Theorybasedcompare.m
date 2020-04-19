   
    D=4*pr.*(1-pr);
    G=exp(gamaset/n);
    ther_err_bsc=D.^3./(1-2*D).^2;
    ther_err_gau=G.^6./(1-(2.*(G).^2)).^2;
    semilogy(gamaset,p_err_bsc,'o-','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
    hold on;
    semilogy(gamaset,p_err_gau,'^-','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
    hold on;
    semilogy(gamaset,p_err_bec,'d-','linewidth',2,'color',[0 0.4 0.9],'markerfacecolor',[0 0.4 0.9],'markeredgecolor',[0 0.7 0.9]);
    hold on;
    semilogy(gamaset,ther_err_bsc,':','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
    hold on;
    semilogy(gamaset,ther_err_gau,'*:','linewidth',2,'color',[0 0.6 0.6],'markerfacecolor',[1 0.5 0.6],'markeredgecolor',[0 0.5 0]);
    xlabel('SNR per Bit in dB'); ylabel('Probability of Bit Error For k=3'); grid
    legend('BSC','Gaussian Noise','BEC','ther-BSC','ther-GAU');     
