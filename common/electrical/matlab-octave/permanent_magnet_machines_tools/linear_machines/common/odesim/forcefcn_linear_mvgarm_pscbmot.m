function [FTVec, FA, ForceBD, xRmachine, xRmagcouple, vR, unitv] = forcefcn_linear_mvgarm_pscbmot(design, simoptions, xT, vT, xA, vA, xBh, vBh, xBs, vBs)

    xR = ((xT + simoptions.xEoffset) - xA);
    
    xRmachine = xR / design.PoleWidth;
    
    xRmagcouple = xR / design.MagCouple.Wr;
    
    vR = vT - vA;
    
    FEAFy = evalhpslmodd(design.MagCouple.slm_FEAFy, 0, xRmagcouple, 1) ...
                * design.MagCouple.N * simoptions.NoOfMachines;
    
%     FEAFy = evalhpslmodd(design.slm_FEAFy, 0, xR, 1) * design.Poles(1) *
%     design.sides * simoptions.NoOfMachines;
    
    % Add a drag force with linear relationship to the relative velocity
    % (useful for simulating eddy current forces etc.). These will be
    % calculated as positive for positive relative velocity of the
    % translator relative to the armature, i.e. vR = vT - vA.
    FLinearDrag = design.klineardrag * vR; % * abs(cos(pi*xR)); % design.klineardrag * (vT - vA); %
    
    % calculate the translator friction
    FfT = friction(design.mu_fT, 9.81 * design.massT .* sin(design.AngleFromHorizontal)) * -sign(vT);
    
    % calculate the translator drag
    FdragT = -sign(vT) .* 0.5 .* realpow(vT,2) .* simoptions.DragRho ...
                .* design.Cd .* design.EffDragArea .* design.NStages .* simoptions.NoOfMachines;
    
    FT = FfT + FdragT - FEAFy - FLinearDrag;
    
    % Find unit vector in the direction pointing from hawse hole to the buoy
    unitv = [simoptions.BuoySim.tether_length+xBh, xBs] / norm([simoptions.BuoySim.tether_length+xBh, xBs]);
    
    % Get the force exerted on the buoy
    FTVec = FT * unitv;
    
    % determine the spring forces on the armature at position xA
    Fs = Fs_Snapper(xA - design.xSC, design);
        
    % Now determine the frictional forces on the armature

%     % The weight of the armature resolved according to the angle form being
%     % mounted vertically
%     N = 9.81 * design.massA .* sin(design.AngleFromHorizontal);
%     
%     % Calculate the magnitude of the frictional force when moving
%     FfA = friction(design.mu_fA, N);
    
    % Calculate the drag forces on the armature
    FdragA = -sign(vA) .* 0.5 .* realpow(vA,2) .* simoptions.DragRho ...
                .* design.Cd .* design.ReacDragArea .* design.NStages .* simoptions.NoOfMachines;

    % return the forces on the reactor for the purposes of calculating its
    % acceleration
    FA = [Fs, FdragA, FEAFy, FLinearDrag];
    
    % force breakdown, of all the machine forces, contains:
    % [Fs, FdragA, FEAFy, FLinearDrag, FfT, FdragT]
    ForceBD = [FA, FfT, FdragT];
 
end