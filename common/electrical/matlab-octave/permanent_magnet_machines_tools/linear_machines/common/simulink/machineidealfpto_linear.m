function [EMF, idealIcoils, limitIcoils, Force, Ploss] = machineidealfpto_linear(design, simoptions, Fpto, xBh, xBs, vBh, vBs, Icoils)

    % first determine xT from the new tether length, change in translator
    % vertical position will be change in distance from hawser to buoy, i.e
    % the change in tether length
%     xT = sqrt((xBh + simoptions.BuoySim.tether_length).^2 + xBs.^2) - simoptions.BuoySim.tether_length;
    xT = xBh;

    % convert to position relative to pole width
    xR = xT ./ design.PoleWidth;
    
    % Get the positions of the three coils in a 3-phase block based on the
	% translator position
    pos = -[xR-(2/3); xR; xR+(2/3)];
    
    % Find dpsidxR from an slm object
    dpsidxR = slmpsidot_slotless(design, pos, design.PoleWidth);
    
%     % Find unit vector in the direction pointing from hawse hole to the buoy
%     unitv = [simoptions.BuoySim.tether_length+xBh, xBs] / norm([simoptions.BuoySim.tether_length+xBh, xBs]);
% 
%     % Then find dot product of heave and surge velocities with unit vector
%     % to get correct direction and magnitude of translator velocity
%     vT = dot(unitv, [vBh, vBs]);

    vT = vBh;

    % determine the emf (voltage) in the coils, - vR * dpsi / dx
    EMF = - vT .* dpsidxR;
    
    % determine what currents are required to obtain the desired power
    % take-off force by vector division of the desired PTO force by the
    % rate of change in flux linkage w.r.t. x
    % The force on the armature is the opposite to the force on the
    % translator ( the power take-off force )
    root2 = sqrt(2);
    
%     Irms = -(Fpto / (design.Poles(1) * simoptions.NoOfMachines)) .* root2 / (3 .* slmpsidot_slotless(design, 0.5, design.PoleWidth));
    
    Irms = -(Fpto / (design.Poles(1) * simoptions.NoOfMachines)) .* root2 / design.Maxdlambdadx;
    
%     Icoils = (Irms .* root2) * (dpsidxR / norm(dpsidxR)); %  put in for the old method

    idealIcoils = (Irms .* root2) * (dpsidxR / (-design.Maxdlambdadx/3));

    Force = -sum(idealIcoils .* -dpsidxR) .* design.Poles(1) .* simoptions.NoOfMachines;

    scale = Fpto / Force;
    
    idealIcoils = idealIcoils .* scale;
    
%     idealIcoils = vT * Fpto * unit(dpsidxR) / dpsidxR;
    
    Irms = Irms * scale;
    
    Jrms = Irms / design.ConductorArea;
    
    if abs(Jrms) > simoptions.maxAllowedJrms
        
        Irms = sign(Irms) * simoptions.maxAllowedJrms .* design.ConductorArea;
        
%         limitIcoils = (Irms .* root2) * (dpsidxR / norm(dpsidxR)); %  put in for the old method

        limitIcoils = (Irms .* root2) * (dpsidxR / (-design.Maxdlambdadx/3));
        
    else
        
        limitIcoils = Icoils;
        
    end
    
    % determine the forces due to the magnets and electrical forces at
    % the relative position xR absed on the coil current and rate of change
    % of flux linkage w.r.t. xR at this point.
    Force = -sum(Icoils .* -dpsidxR) .* design.Poles(1) .* simoptions.NoOfMachines;
    
    %Ploss = sum(Icoils.^2 .* design.CoilResistance) .* design.Poles(1) .* simoptions.NoOfMachines;
    Ploss = sum(Icoils.^2 .* design.CoilResistance) .* design.Poles(1) .* simoptions.NoOfMachines;

end