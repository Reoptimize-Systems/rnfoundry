function [EMF, torque, flux, design] = rotgensfcn_common(design, simoptions, theta, omega, phaseI)
% Calculates the EMF and torque from a single rotary generator given the
% phase currents, rotational speed and position.
%
% Syntax
%
%   [EMF, torque] = rotgensfcn_common(design, simparams, displacement,...
%                                       velocity, phaseI)
% 
% Input
% 
%   design       - design structure
%   simparams    - simparams structure
%   theta        - angular displacement in rad
%   omega        - angular velocity in rad/s
%   current      - a vector of pahse currents
% 
% Output
% 
%   EMF     - a vector of phase EMFs
%   Torque  - a scalar of the torque due to the coil currents and losses
% 
% See also: rotgen_three_phase_sfcn
 
% Created Richard Crozier 2013

    % get the coil currents from the phase currents
    Icoils = phaseI ./ design.Branches;
    
    % Calculate the EMF and force at the reference point
    [Tqeff, Tqreac, EMF, dpsidthetaR, design, pos] = ...
        machineodesim_AM(design, simoptions, theta, 0, omega, 0, Icoils);
    
    EMF = EMF .* design.CoilsPerBranch;
    
    % get the additional torques due to losses
    TqaddE = losstorques_ROTARY(design, simoptions, theta, omega);
    
    torque = (Tqeff + TqaddE);
    
    % total phase flux linkage
    flux = periodicslmeval(pos, design.slm_fluxlinkage, 0, false) .* design.CoilsPerBranch .* design.Branches;
    
end