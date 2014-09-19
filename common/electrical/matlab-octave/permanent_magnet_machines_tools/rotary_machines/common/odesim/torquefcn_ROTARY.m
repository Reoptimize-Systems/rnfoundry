function [Torque, TorqueBD] = torquefcn_ROTARY(design, simoptions, thetaE, omegaE, EMF, Iphases)
% calculates forces for a linearised rotary machine dynamic simulation
%
% Syntax
%
% [Torque, TorqueBD] = torquefcn_ROTARY(design, simoptions, thetaE, omegaE, EMF, Iphases)
%
% 

    thetaR = thetaE ./ design.PoleWidth;
    
    [~, TqLiron, TqLeddy] = losstorques_ROTARY(design, simoptions, thetaR, omegaE);
    
    TorqueBD = [0, 0, TqLiron, TqLeddy];
    
    % sum the forces
    Torque = sum(TorqueBD);

end