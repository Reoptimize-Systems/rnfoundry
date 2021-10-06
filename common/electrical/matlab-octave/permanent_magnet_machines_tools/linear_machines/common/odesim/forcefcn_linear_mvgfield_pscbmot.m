function [FT, FF] = forcefcn_linear_mvgfield_pscbmot(design, simoptions, xT, vT, xF, vF, xBh, xBs, vBh, vBs)

    xR = (xT - xF);
    
    xRmachine = xR / design.PoleWidth;
    
    xRmagcouple = xR / design.MagCouple.Wr;
    
    FEAFy = evalhpslmodd(design.MagCouple.slm_FEAFy, 0, xRmagcouple, 1) * design.MagCouple.N * simoptions.NoOfMachines;
    
    % Add a drag force with linear relationship to the relative velocity
    % (useful for simulating eddy current forces etc.). These will be
    % calculated as positive for positive relative velocity of the
    % translator relative to the armature, i.e. vR = vT - va.
    FLinearDrag = design.klineardrag * (vT - vF) * abs(sin(pi*xRmachine)); % design.klineardrag * (vT - vA); %
    
    % calculate the translator friction
    FfT = friction(design.mu_fT, 9.81 * design.massT .* sin(design.AngleFromHorizontal)) * -sign(vT);
    
    FT = FEAFy + FLinearDrag - FfT;
    
    % Find unit vector in the direction pointing from hawse hole to the buoy
    unitv = [simoptions.BuoySim.tether_length+xBh, xBs] / norm([simoptions.BuoySim.tether_length+xBh, xBs]);
    
    % Get the force exerted on the buoy
    FT = -FT * unitv;
    
    % determine the spring forces on the field at position xF
    Fs = Fs_Snapper(xF - design.xSC, design);
        
    % Now determine the frictional forces on the field

%     % The weight of the field resolved according to the angle form being
%     % mounted vertically
%     N = 9.81 * design.massF .* sin(design.AngleFromHorizontal);
%     
%     % Calculate the magnitude of the frictional force when moving
%     Ffricfield = friction(design.mu_fF, N);
%     Ffricfield = 0;
    
    % Calculate the drag forces on the field
    FdragF = -sign(vF) .* 0.5 .* realpow(vF,2) .* simoptions.BuoySim.SeaParameters.rho .* design.Cd .* design.DragArea;

    % total the force on the reactor
    FF = [Fs, FLinearDrag, FdragF, FEAFy];
 
end