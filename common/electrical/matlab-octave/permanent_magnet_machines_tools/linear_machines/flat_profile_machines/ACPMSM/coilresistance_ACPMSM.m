function R = coilresistance_ACPMSM(design)
% CoilResistance_Snapper: calculates the resistance of a single coil
% winding based on the secifications in the design structure
%

    % there is assumed to be an arc at the corners with radius half the
    % coil cross-section width
    length = design.CoilTurns .* rectcoilmtl(design.ls, design.Taup - design.Wc, design.Wc);
    
    R = 1.68e-8 * length / (pi * (design.Dc/2)^2);

end