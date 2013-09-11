function [design, simoptions] = finfun_RADIAL(design, simoptions)
% perform post-processing common to all TORUS type machines of the data
% produced by an fea simulation in readiness for a dynamic simulation

    % do post-processing common to TORUS type machines
%     
%     % Convert air gap closing force to force per unit area
%     design.gforce = design.gforce ./ (design.taupm * design.hm);
%     
%     % fit a polynomial to the air gap force
%     design.p_gforce = polyfitn(design.gvar, design.gforce, 2);

    % set the polewidth to be the distance in radians swept out by a pole
    design.PoleWidth = design.thetap;
    
    % call finfun_ROTARY
    [design, simoptions] = finfun_ROTARY(design, simoptions);

end 