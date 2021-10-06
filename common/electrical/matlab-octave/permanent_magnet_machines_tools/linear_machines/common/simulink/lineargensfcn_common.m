function [EMF, force, flux, design] = lineargensfcn_common(design, simoptions, xR, vR, phaseI)
% Calculates the EMF and torque from a single rotary generator given the
% phase currents, rotational speed and position.
%
% Syntax
%
%   [EMF, force] = lineargensfcn_common (design, simparams, displacement,...
%                                        velocity, phaseI)
% 
% Input
% 
%   design       - design structure
%   simparams    - simparams structure
%   xR           - relative displacement in m of stator and translator
%   vR           - relative velocity in m/s of stator and translator
%   current      - a vector of pahse currents
% 
% Output
% 
%   EMF     - a vector of phase EMFs
%   force   - a scalar of the torque due to the coil currents and losses
% 
% See also: rotgen_three_phase_sfcn
 
% Created Richard Crozier 2013

    % get the coil currents from the phase currents
    Icoils = phaseI ./ design.Branches;
    
    % Calculate the EMF and force at the reference point
    [Feff, ~, EMF, ~, design, pos] = ...
        machineodesim_AM (design, simoptions, xR, 0, vR, 0, Icoils);
    
    EMF = EMF .* design.CoilsPerBranch;
    
    % get the additional torques due to losses
    FaddE = lossforces_AM (design, simoptions, xR, vR);
    
    force = (Feff + FaddE);
    
    % total phase flux linkage
    flux = periodicslmeval (pos, design.slm_fluxlinkage, 0, false) .* design.CoilsPerBranch .* design.Branches;
    
end