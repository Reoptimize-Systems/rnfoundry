function [dpsidxRF, EMF, FF, ForceVec, xE, vE, design] = machineodesim_linear_mvgfield(design, simoptions, Icoils, xF, vF, xBh, xBs, vBh, vBs)
% machineodesim_linear_mvgfield: determines the EMF, rate of change in flux
% linkage with relative armature/translator position, Shear Force and the
% velocities of each part of the machine

    %**********************************************************************
    % Determines the forces acting within the machine, and other important
    % machine outputs

    % next we determine the relative position of the armature and field in
    % order to determine forces and flux linkage

    % first determine xT from the new tether length, change in translator
    % vertical position will be change in distance from hawser to buoy, i.e
    % the change in tether length
    xE = (sqrt((xBh + simoptions.BuoySim.tether_length).^2 + xBs.^2) - simoptions.BuoySim.tether_length);

    % Find unit vector in the direction pointing from hawse hole to the buoy
    unitv = [simoptions.BuoySim.tether_length+xBh, xBs] / norm([simoptions.BuoySim.tether_length+xBh, xBs]);

    % Then find dot product of heave and surge velocities with unit vector
    % to get correct direction and magnitude of translator velocity
    vE = dot(unitv, [vBh, vBs]);

    % Get the EMF and forces using the core machine simulation function
    [FF, FA, EMF, dpsidxRF, design] = machineodesim_AM(design, simoptions, xF, xE, vF, vE, Icoils);

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
    ForceVec = FA * unitv;
    
end