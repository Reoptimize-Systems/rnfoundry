function [Force, ForceBD, xR] = forcefcn_linear_pscbmot(design, simoptions, xT, vT, EMF, Iphases)

    % calculate the displacement relative to the pole width
    xR = xT ./ design.PoleWidth;
    
    % Add a drag force with linear relationship to the relative velocity
    % (useful for simulating eddy current forces etc.). These will be
    % calculated as positive for positive relative velocity of the
    % translator relative to the armature, i.e. vR = vT - va.
    FLinearDrag = design.klineardrag * -vT; %
    
    % calculate the translator friction
    FfT = friction(design.mu_fT, 9.81 * design.massT .* sin(design.AngleFromHorizontal)) * -sign(vT);
    
    % calculate the forces due to losses
    FLoss = lossforces_AM(design, simoptions, xR, vT);

    ForceBD = [FLinearDrag, FfT, FLoss];
    
    % return the force exerted on the prime mover
    Force = sum(ForceBD);
    
end