function R0 = tR0_PMSM(qvc, rho_c, kthc)
    
    %Equiv thermal resistance of cooling duct

    R0= 1/(qvc.*rho_c.*kthc);


end