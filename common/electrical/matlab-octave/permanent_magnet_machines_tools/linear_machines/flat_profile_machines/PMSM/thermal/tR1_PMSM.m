function R1 = tR1_PMSM(ls, bt, alpha_1)

    % Equiv thermal resistance of stator back (tooth width) convection to
    % cooling air
    R1 = (3.*ls*bt.*alpha_1).^-1; 
    
end