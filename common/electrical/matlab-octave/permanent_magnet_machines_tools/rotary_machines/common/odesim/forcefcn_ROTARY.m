function [Force, ForceBD] = forcefcn_ROTARY(design, simoptions, xE, vE, EMF, Iphases)
% calculates forces for a linearised rotary machine dynamic simulation
%
% Syntax
%
% Force = forcefcn_ROTARY(design, simoptions, xT, vT, xBh, xBs, vBh, vBs)
%
% 

    xR = xE ./ design.PoleWidth;
    
    [Flosstot, Flossiron, Flosseddy] = lossforces_ROTARY(design, simoptions, xR, vE);
    
    ForceBD = [0, 0, Flossiron, Flosseddy];
    
    % sum the forces
    Force = sum(ForceBD);

end