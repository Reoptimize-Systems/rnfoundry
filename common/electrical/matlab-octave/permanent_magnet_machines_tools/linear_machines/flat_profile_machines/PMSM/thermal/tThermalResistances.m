function R = tThermalResistances(q_vcair, rho_cair, k_thcair, k_Cu, bt, bs, hsy, hry, ls, dsi, hm0, hm1, hm, hi, hs, hs1, hs2, Taup, bs1, bm, Nphases, pp, Tautb, Q, alpha, lambda)

    bCu = bs-2 * hi;

    hCu = (hs / 2) - 4 * hi;

    lb = 2 * Taup; %end winding
    
    Taus = Taup / Nphases;

    %Thermal resistances from Grauers Thesis
    
    %Equiv thermal resistance of cooling duct, R_0
    R(1) = (q_vcair * rho_cair * k_thcair)^-1;
    %Equiv thermal resistance of stator back (tooth width) convection to
    %cooling air, R_1
    R(2) = (3 * ls * bt * alpha(1))^-1; 
    %Equiv thermal resistance of stator back (slot width) convection to
    %cooling air, R_2
    R(3) = (3 * ls * bs * alpha(1))^-1; 

    % In Grauers R_3, R_4, R_5, R_6, R_12, R_14, R_21 use l_u not ls -
    % difference l_u = ls * k_Fe where k_Fe is iron fill factor
    
    %Equiv thermal resistance of half stator back (tooth width) yoke radial
    %conduction
    R(4) = (0.5 * hsy) / (ls * bt * lambda(1)); 
    %Equiv thermal resistance of half stator back (slot width) yoke radial
    %conduction
    R(5) = (0.5 * hsy) / (ls * bs * lambda(1)); 
    %Equiv thermal resistance of half tooth yoke tangential conduction
    R(6) = (0.5 * bt) / (ls * hsy * lambda(1)); 
    %Equiv thermal resistance of half slot yoke tangential conduction
    R(7) = (0.5 * bs) / (ls * hsy * lambda(1)); 
    %mean adjustment for R(5 heat
    R(8) = R(6) / -3; 
    %mean adjustment for R_6 heat
    R(9) = R(7) / -3; 
    %mean adjustment for R_3 hea
    R(10) = R(4) / -3; 
    %mean adjustment for R_4 heat
    R(11) = R(5) / -3; 
    
    % Equiv thermal resistance of insulation (radial from slot) conduction
    % bCu is conductor width bCu = bs-2 * hi
    R(12) = hi / (ls * bCu * lambda(4)); 
    % Equiv thermal resistance of tooth (radial from level with top
    % conductor) conduction %hCu is conductor height hCu = (hs / 2)-4 * hi
    R(13) = 0.5 * (hCu + 2 * hi) / (ls * bt * lambda(1)); 
    % Equiv thermal resistance of half conductor (radial) conduction
    R(14) = (0.5 * hCu) / (ls * bCu * lambda(2)); 
    % Equiv thermal resistance of half tooth (centre tooth to edge of slot)
    % conduction
    R(15) = (0.5 * bt) / (ls * (hCu + 2 * hi) * lambda(1)); 
    % Equiv thermal resistance of insulation (radial from conductor)
    % conduction
    R(16) = hi / (ls * hCu * lambda(4)); 
    % Equiv thermal resistance of conductor (tangential) conduction
    R(17) = (0.5 * bCu) / (ls * hCu * lambda(2)); 
    % mean adjustment for P_c R_14 heat
    R(18) = R(15) / -3; 
    % mean adjustment for P_d R_16 heat
    R(19) = R(17) / -3; 
    % mean adjustment for P_c R_12 heat
    R(20) = R(13) / -3; 
    % mean adjustment for P_d R_13 heat
    R(21) = R(14) / -3; 
    
    % Equiv thermal resistance radial tooth widening conduction
    R(22) = (hs1 + hs2) / (ls * (0.5 * (Tautb + bt)) * lambda(1)); 
    % Equiv thermal resistance radial tooth to airgap convection
    R(23) = 1 / (ls * (Taus-bs1) * alpha(2)); 
    % Equiv thermal resistance radial airgap to magnet convection  +
    % radial conduction GRP band
    R(24) = (1 / (ls * bm * alpha(2))) + (hm1 / (ls * bm * lambda(7))); 
    % Equiv thermal resistance radial half PM conduction
    R(25) = (0.5 * hm) / (ls * bm * lambda(5)); 
    % mean adjustment for P_e R_24 heat
    R(26) = R(25) / -3; 
    % Equiv thermal resistance radial PM glue conduction
    R(27) = hm0 / (ls * bm * lambda(6)); 
    % Equiv thermal resistance radial rotor yoke conduction
    R(28) = hry / (ls * Taup * lambda(1)); 
    % Equiv thermal resistance radial rotor yoke convection to air
    R(29) = 1 / (ls * Taup * alpha(5)); 
    
    % Change for linear: 
    
    %Equiv thermal resistance end shield convection to air
    R(30) = 1 / (pi * ((0.5* dsi) + hs + hsy)^2 * alpha(3));
    % Equiv thermal resistance of half axial conductor length conduction
    R(31) = 0.5 * ls / (hCu * bCu * k_Cu * lambda(3)); 
    
    R(32) = R(31);
    %mean adjustment for P_d R_31 heat
    R(31) = R(32) / -3; 
    %Equiv thermal resistance of half end winding length conduction %lb is
    %end winding length, lb = 2 * Taup
    R(33) = 0.5 * lb / (hCu * bCu * k_Cu * lambda(3)); 
    %mean adjustment for P_f R_32 heat
    R(34) = R(33) / -3; 
    %Equiv thermal resistance of half end winding radial conduction
    R(37) = 0.5 * hCu / (lb * bCu * lambda(2)); 
     %mean adjustment for P_f R_36 heat
    R(35) = R(37) / -3;
     %Equiv thermal resistance of end winding radial convection
    R(36) = 2 / (lb * bCu * alpha(4));
    %Equiv thermal resistance of end winding axial convection
    R(39) = 2 / (lb * hCu * alpha(4)); 
    %Equiv thermal resistance of half end winding axial conduction
    R(40) = 0.5 * bCu / (lb * hCu * lambda(2)); 
    %mean adjustment for P_f R_39 heat
    R(38) = R(40) / -3; 

    % Simplified model thermal resistances from Grauers
    R(51) = R(2);
    R(52) = (R(3) + R(4)) / Q;
    R(53) = (R(4) + R(5)) / Q;
    R(54) = (R(7) + R(9) + R(10) + R(11) + 0.5 * (R(6) + R(7))) / Q;
    R(55) = (R(4) + R(13)) / Q;
    R(56) = (R(5) + R(12) + R(14)) / Q;
    R(57) = ((R(20) + R(18) + R(19)) + 0.5 * (R(15) + R(16) + R(17))) / Q;
    R(58) = R(21) / Q;
    R(59) = (2 * R(13)) / Q;
    R(60) = (2 * R(14) + 2 * R(12) + R(21)) / Q;

    R(61) = ((R(13) + R(22) + R(23)) / Q) + ((R(24) + R(25)) / (2 * pp));
    R(62) = (R(25) + R(27) + R(28) + R(29)) / 2 / pp;
    R(63) = R(28);
    R(64) = (R(31) + 0.5 * (R(32) + R(33))) / Q;

    R_64a = 0.5 * (R(35) + (0.5 * (R(37) + R(36)))) / Q;
    R_64b = 0.5 * (R(36) + (0.5 * (R(40) + R(39)))) / Q;

    R(65) = (R_64a^-1 + R_64b^-1)^-1;
    R(66) = 0.5 * R(34) / Q;

end