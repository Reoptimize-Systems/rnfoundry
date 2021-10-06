function Pmat = tNodalPowerMatrix_PMSM(Tamb, P_Hyt, P_Hysy, P_Ect, P_Ecsy, P_Ecm, P_ad, P_Cus, ls, lb, bt, bs, Taus, G_50, G_62)
   
    %P_n  =  Heat input at node n
    
    P_0 = 0; 
    
    P_1 = (bt / Taus) * (P_Hysy + P_Ecsy);
    
    P_2 = (bs / Taus) * (P_Hysy + P_Ecsy);
    
    P_3 = 0.5 * (P_Hyt + P_Ect);
    
    P_5 = 0;
    
    P_6 = 0;
    
    P_8 = 0;
    
    P_9 = P_3 + P_ad;
    
    P_11 = P_Ecm;
    
    P_4 = 0.5 * (ls / (ls + lb)) * P_Cus;
    
    P_7 = 0.5 * (lb / (ls + lb)) * P_Cus;
    
    P_10 = P_4;
    
    %Power matrix, [P_mat] 1by12
    Pmat = [(P_0 + (Tamb * G_50)); P_1; P_2; P_3; P_4; P_5; P_6; P_7; (P_8 + (Tamb * G_62)); P_9; P_10; P_11];

end