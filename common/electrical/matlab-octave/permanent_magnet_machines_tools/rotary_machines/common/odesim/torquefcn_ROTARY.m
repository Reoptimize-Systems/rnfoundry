function [Torque, TorqueBD] = torquefcn_ROTARY(design, simoptions, thetaE, omegaE, EMF, Iphases)
% calculates forces for a linearised rotary machine dynamic simulation
%
% Syntax
%
% Force = forcefcn_ROTARY(design, simoptions, xT, vT, xBh, xBs, vBh, vBs)
%
% 

    thetaR = thetaE ./ design.PoleWidth;
    
    [TqLtot, TqLiron, TqLeddy] = losstorques_ROTARY(design, simoptions, thetaR, omegaE);
    
    TorqueBD = [0, 0, TqLiron, TqLeddy];
    
    % sum the forces
    Torque = sum(TorqueBD);

end