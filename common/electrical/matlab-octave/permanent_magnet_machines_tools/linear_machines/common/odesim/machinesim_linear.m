function [dpsidxF, EMF, Feff, ForceVec, xE, vE, unitv] = machinesim_linear(design, simoptions, Icoils, xBh, xBs, vBh, vBs)
% determines the EMF, rate of change in flux linkage with relative
% armature/translator position, Shear Force and the velocities of each part
% of the machine
%
% Syntax
%
% [dpsidxF, EMF, Feff, ForceVec, xE, vE] = ...
%       machinesim_linear(design, simoptions, Icoils, xBh, xBs, vBh, vBs)
%
% Input
%
% design and simoptions are structures which include certain fields
% necessary for the machine simulation. design must contain at least the
% following fields:
%   
%   FieldDirection - a scalar value of either 1 or -1 indicating the
%       direction of the magnetic field displacement relative to the
%       direction of the prime mover displacement. If 1 this means the
%       magnetic field moves in the same direction as the translator, i.e.
%       the magnets are mounted on the translator and the coils are
%       stationary. If -1, the opposite is true, and the magnetic field is
%       moving in the opposite direction.
%
%   PoleWidth - scalar value of length of a pole, or half period of the
%       flux waveform
%
%   slm_psidot - An slm object fitted to the flux linkage from the point of
%       maximum flux linkage to the point of minimum flux linkage over one
%       half period of the flux linkage waveform. The displacement over the
%       half period should be normalised to the domain 0 to 1.
%
%   phases - scalar number of phases in the machine armature.
%
%   PowerPoles - scalar value of the number active poles in the machine
%       armature.
%
% simoptions must contain at least the following fields:
%
%   tether_length - this is the starting lenght of a rigid tether from the
%       prime mover to the effector of the machine.
%
%   maxAllowedxT - scalar value of the maximum allowed excursion of the
%       effector from its initial position. This is used to apply end
%       stops. Can be set to inf, but must be present.
%
%   EndStopks - this is the spring constant of the end stop encountered at
%       the maxAllowedxT position
%
%   NoOfMachines - the number of machines linked together for the purposes
%       of force calculation.
%
% The following arguments are then:
%
% Icoils - A vector of values of the the coil current, should be the same
%   length as the number of phases in the machine.
% 
% xBh - The displacement of the prime mover (typically a buoy) in the
%   vertical (heave) direction.
%
% xBs - The displacement of the prime mover (typically a buoy) in the
%   horizontal (surge) direction.
% 
% vBh - The velocity of the prime mover (typically a buoy) in the
%   vertical (heave) direction.
% 
% vBs - The velocity of the prime mover (typically a buoy) in the
%   horizontal (surge) direction.
%
% Output
% 
% dpsidxF -  vector of values of the derivative of the flux linkage w.r.t.
%   displacement at the current position of the machine.
% 
% EMF - vector of values of the EMF in a single coil in each phase
% 
% Feff - The force acting on the effector of the machine (the part attached
%   to the prime mover). This force is defined as positive in the direction
%   of positive displacement of the effector.
% 
% ForceVec - A vector of value of the force acting on the prime mover in
%   the surge and heave directions
% 
% xE -  The vertical displacement of the effector of the machine.
% 
% vE - The vertical velocity of the effector of the machine.
%


    % first determine the translator displacement from the tether
    % length, the change in translator vertical position will be the change
    % in distance from hawser to buoy, i.e the change in tether length
    xE = sqrt((xBh + simoptions.tether_length).^2 + xBs.^2) - simoptions.tether_length;

    xE = xE + simoptions.xEoffset;
    
    % Find unit vector in the direction pointing from hawse hole to the
    % buoy, we add a tiny length to the tether length in case it is of
    % length zero which would make calculation of the unit vector
    % impossible
    unitv = [simoptions.tether_length+1e-6+xBh, xBs] / norm([simoptions.tether_length+1e-6+xBh, xBs]);

    % Then find dot product of heave and surge velocities with unit vector
    % to get correct direction and magnitude of the velocity of the
    % magnetic field relative to the coils
    vE = dot(unitv, [vBh, vBs]);

    % Get the EMF and forces using the core machine simulation function
    [Feff, Freac, EMF, dpsidxF] = coremachinesim_linear(design, simoptions, xE, 0, vE, 0, Icoils);

    % If the translator is at the maximum allowed displacement and
    % continuing to move in that direction, add a retarding force to slow
    % it down, which is related to the amount by which the max allowed
    % displacement has been exceeded
%     if (xT > simoptions.maxAllowedxT && vT > -simoptions.EndStopVThresh) || (xT < -simoptions.maxAllowedxT && vT < simoptions.EndStopVThresh)
% 
%         Feff = Feff - sign(xT) * simoptions.EndStopks * (abs(xT) - simoptions.maxAllowedxT);
%         
%     end
    
    
    %  Calculating the proportion that the machine forces are producing on
    %  surge and in heave. Do this by multiplying unit vector in direction
    %  of hawser to buoy by total forces. This results in new vector:
    %
    %               FfeaVec = [Ffea_heave, Ffea_surge]
    %
    %  Which should also have forces in the right directions. The FEA
    %  forces are calculated for the coils, therefore when the translator
    %  is moving up, Ffea is also acting upwards. We must reverse their
    %  direction to get the forces acting on the translator (and buoy).
    ForceVec = Feff * unitv;
    
end