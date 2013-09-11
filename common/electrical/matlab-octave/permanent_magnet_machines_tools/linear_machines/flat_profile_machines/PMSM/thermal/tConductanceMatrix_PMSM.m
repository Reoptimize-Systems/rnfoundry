function Gmat = tConductanceMatrix_PMSM(Rs)

    % invert the thermal resistances to get conductances
    Rs = Rs .^ -1;
    %Thermal Conductances

    G_51 = Rs(1);
    G_52 = Rs(2);
    G_53 = Rs(3);
    G_54 = Rs(4);
    G_55 = Rs(5);
    G_56 = Rs(6);
    G_57 = Rs(7);
    G_58 = Rs(8);
    G_59 = Rs(9);
    G_60 = Rs(10);
    G_61 = Rs(11);
    G_62 = Rs(12);
    G_63 = Rs(13);
    G_64 = Rs(14);
    G_65 = Rs(15);
    
    %Conductance matrix, [G_mat], 12by12
    G_00 = G_51 + G_52 + G_53;
    G_11 = G_51 + G_53 + G_54;
    G_22 = G_52 + G_53 + G_55;
    G_33 = G_54 + G_56 + G_58;
    G_44 = G_56 + G_57 + G_63;
    G_fivefive = G_55 + G_57 + G_59;
    G_sixsix = G_63 + G_65 + G_63;
    G_77 = G_65 + G_64;
    G_88 = G_62 + G_64 + G_61;
    G_99 = G_58 + G_56 + G_60;
    G_1010 = G_56 + G_59 + G_63;
    G_1111 = G_60 + G_61;

    Gmat = [G_00 -1 * G_51 -1 * G_52 0 0 0 0 0 0 0 0 0;
            -1 * G_51 G_11 -1 * G_53 -1 * G_54 0 0 0 0 0 0 0 0;
            -1 * G_52 -1 * G_53 G_22 0 0 -1 * G_55 0 0 0 0 0 0;
            0 -1 * G_54 0 G_33 -1 * G_56 0 0 0 0 -1 * G_58 0 0;
            0 0 0 -1 * G_56 G_44 -1 * G_57 -1 * G_63 0 0 0 0 0;
            0 0 -1 * G_55 0 -1 * G_57 G_fivefive 0 0 0 0 -1 * G_59 0;
            0 0 0 0 -1 * G_63 0 G_sixsix -1 * G_65 0 0 -1 * G_63 0;
            0 0 0 0 0 0 -1 * G_65 G_77 -1 * G_64 0 0 0;
            0 0 0 0 0 0 0 -1 * G_64 G_88 0 0 -1 * G_61;
            0 0 0 -1 * G_58 0 0 0 0 0 G_99 -1 * G_56 -1 * G_60;
            0 0 0 0 0 -1 * G_59 -1 * G_63 0 0 -1 * G_56 G_1010 0;
            0 0 0 0 0 0 0 0 -1 * G_61 -1 * G_60 0 G_1111];


end