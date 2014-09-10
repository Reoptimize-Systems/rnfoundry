function [Torque, TorqueBD] = torquefcn_RADIAL(design, simoptions, thetaE, omegaE, EMF, Iphases)
% calculates forces for a linearised rotary machine dynamic simulation
%
% Syntax
%
% [Torque, TorqueBD] = torquefcn_RADIAL(design, simoptions, thetaE, omegaE, EMF, Iphases)
%
% 

    [Torque, TorqueBD] = torquefcn_ROTARY(design, simoptions, thetaE, omegaE, EMF, Iphases);

end