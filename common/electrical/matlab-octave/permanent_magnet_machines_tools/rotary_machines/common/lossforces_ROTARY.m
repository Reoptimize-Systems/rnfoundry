function [FLtot, FLiron, FLeddy] = lossforces_ROTARY(design, simoptions, xR, vR)
% calculates the forces due to eddy current losses in a rotary machine at a
% given velocity

    [FLtot, FLiron, FLeddy] = lossforces_AM(design, simoptions, xR, vR);
    
%     if vT 
%         
%         % calculate the frequency of oscillation at the core surface m/s / m
%         freq = abs(vT ./ (2 * design.taupm));
% 
%         % determine the skin depth
%         d = skindepth(design.CoreResistivity, design.CoreMuR, freq);
% 
%         % calculate the core losses from empirical data
%         Ploss = d * corelossarea * design.CoreMaterialDensity ...
%                 * interp2(simoptions.CoreLossData.fq, ...
%                           simoptions.CoreLossData.Bq, ...
%                           simoptions.CoreLossData.Pq, ...
%                           freq, ...
%                           design.MaxCoreSurfaceField, ...
%                           '*cubic');
%         
%         % force on one side of the core is then power divided by velocity
%         FaddE = -Ploss / vT;
%         
%         % total force is both sides times the number of stages
%         FaddE = FaddE * 2 * design.NStages;
%     
%     else
%         
%         % if velocity is zero, force is zero
%         FaddE = 0;
%         
%     end
    
end