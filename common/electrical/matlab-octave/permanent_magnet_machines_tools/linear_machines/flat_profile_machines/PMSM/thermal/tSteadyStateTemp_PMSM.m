function [Tmat, Rs] = tSteadyStateTemp_PMSM(Gmat, Tamb, Is, R, rho_Cu_25, alpha_Cu, P_Hyt, P_Hysy, P_Ect, P_Ecsy, P_Ecm, P_ad, ls, lb, bt, bs, Taus, l_Cus, A_Cus, Nphases)
    
    Tmat(1:size(Gmat,1),1) = Tamb;
    delTmat = ones(size(Tmat));
    i = 0;
    
    while i < 20 && max(delTmat) > 0.1

        %Mean conductor temperature, dependent on I and T
        Tmean = ((2*ls/l_Cus) * 0.5 * (Tmat(5) + Tmat(11))) + ((2*lb/l_Cus) * Tmat(8));

        % The electrical resistance of the conductors depends on temperature
        Rs = l_Cus * rho_Cu_25 * (1 + alpha_Cu * (Tmean - 25)) / A_Cus; %dep on T

        PCus = Nphases * Is^2 * Rs; %dep on T and I_s

        oldTmat = Tmat;

        %Power matrix, [P_mat] 1by12
        %P_mat = [(P_0 + (T_amb * G_50)); P_1; P_2; P_3; P_4; P_5; P_6; P_7; (P_8 + (T_amb * G_62)); P_9; P_10; P_11];
        Pmat = tNodalPowerMatrix_PMSM(Tamb, P_Hyt, P_Hysy, P_Ect, P_Ecsy, P_Ecm, P_ad, PCus, ls, lb, bt, bs, Taus, 1/R(51), 1/R(63));

        % Invert the conductance matrix and multiply by the power matrix to
        % get the temperature matrix, [Tmat], use Gmat \ Pmat, (Matlab
        % matrix left division operation), much faster than doing a straight
        % inversion, but equivalent to inv(Gmat) * Pmat
        Tmat = Gmat \ Pmat; 

        % find the absolute change in temperatures at the last iteration and
        % store in delT_mat
        delTmat = abs(Tmat - oldTmat);

        i = i + 1;

    end
     
end