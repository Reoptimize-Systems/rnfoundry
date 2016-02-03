function [Torque, TorqueBD] = torquefcn_RADIAL_SLOTTED (design, simoptions, thetaE, omegaE, EMF, Iphases)
% calculates forces for a linearised rotary machine dynamic simulation
%
% Syntax
%
% [Torque, TorqueBD] = torquefcn_RADIAL_SLOTTED(design, simoptions, thetaE, omegaE, EMF, Iphases)
%
% 

    [Torque, TorqueBD] = torquefcn_RADIAL (design, simoptions, thetaE, omegaE, EMF, Iphases);
    
    thetaR = thetaE ./ design.PoleWidth;
    
    % calculate the cogging torque
    TorqueBD(end+1) =  periodicslmeval (thetaR, design.slm_coggingtorque, 0, false);
    
    % add it to the total torque
%     Torque = Torque + TorqueBD(end);

end