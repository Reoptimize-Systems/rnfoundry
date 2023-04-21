function [design, simoptions] = finfun_TORUS(design, simoptions)
% perform post-processing common to all TORUS type machines of the data
% produced by an fea simulation in readiness for a dynamic simulation

    % do post-processing common to TORUS type machines
%     
%     % Convert air gap closing force to force per unit area
%     design.ForceGapClosingWithDisp = design.ForceGapClosingWithDisp ./ (design.taupm * design.hm);
%     
%     % fit a polynomial to the air gap force
%     design.p_ForceGapClosingWithDisp = polyfitn(design.DispGapClosingForce, design.ForceGapClosingWithDisp, 2);

    % set the polewidth to be the distance swept out by a pole at the mid
    % point of the magnets
    design.PoleWidth = design.taupm;
    
    % call finfun_ROTARY
    [design, simoptions] = finfun_ROTARY(design, simoptions);

end 