function PMSM_I = TranslatorMomentOfArea_SS_PMSM(IVars)
% TranslatorMomentOfArea_SS_PMSM
%
% A function for calculating the second moment of area of a single pole of
% the single sided linear permanent magnet machine.
%
% Input:
%
%   IVars - (1 x 4) matrix of values for calculating the second moment of
%           area of the teeth and support/central translator section
%
%           IVars(1,1): dbi, half thickness of central section (back 
%                       iron thickness on armature)
%           IVars(1,2): Taup, the pole width
%           IVars(1,3): hs, the tooth height
%           IVars(1,4): bt, the tooth width
%
% Output:
%
%   PMSM_I - second moment of area of single pole for the single sided
%            machine

    hay = IVars(:,1);
    Taup = IVars(:,2);
    hs = IVars(:,3);
    bt = IVars(:,4);

    PMSM_I = ((9 .* (bt.^2) .* (hs.^4)) + ...
              (3 .* Taup .* hs .* bt .* hay .* (hay.^2 + hs.^2)) + ...
              (Taup.^2 .* hay.^2 .* (4 .* hay.^2 + 6 .* hay .* hs + 3 .* hs.^2))) ./ ...
             (12 .* (3 .* bt .* hs + hay .* Taup));
              
%     for i = 1:size(IVars,1)
% 
%         hay = IVars(i,1);
%         Taup = IVars(i,2);
%         hs = IVars(i,3);
%         bt = IVars(i,4);
% 
%         % The global neutral axis will be located at zero in the horizontal
%         % plane, we must therefore first locate the central section axis in
%         % relation to this
%         %Y1 = -0.5*(6*bt*dt + 3*bt*dt*Taup + 3*bt*dt - 2*bc*Taup) / (3*bt*dt + bc*Taup);
% 
%         % Now we will locate the tooth axis in relation to the global neutral
%         % axis
%         %Y2 = Y1 + ((bc+dt)/2);
% 
%         % Calculate second moment of area of the translator central section and
%         % a single tooth
%         I = MomentOfInertiaY1([bc Taup; dt bt], '1.2');
% 
%         % Total second moment of area = I1 + A1Y1^2 + 3(I2 + A2*Y2^2)
%         %PMSM_I = I(1,1) + (Taup*bc*(Y1^2)) + 3 * (I(2,1) + ((Y2^2)*dt*bt)); 
%         
%         % locate the centroid of the section
%         %y1c = translatorY1Centroid(Taup, bc, bt, dt);
%         
%         % calculate the total moment of inertia
%         %PMSM_I(i,1) = I(1,1) + Taup*bc*(y1c - bc/2) + 3.*I(2,1) + 3.*dt.*(y1c - (bc + dt/2));
% 
%         PMSM_I = ((9 .* (bt.^2) .* (hs.^4)) + ...
%                   (3 .* Taup .* hs .* bt .* hay .* (hay.^2 + hs.^2)) + ...
%                   (Taup.^2 .* hay.^2 .* (4 .* hay.^2 + 6 .* hay .* hs + 3 .* hs.^2))) ./ ...
%                   (12 .* (3 .* bt .* hs + hay .* Taup));
%         
%     end  
    
end