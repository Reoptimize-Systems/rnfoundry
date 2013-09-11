function Pfe = ironlosses_PMSM(design, simoptions, vT)
% calculates iron losses in the armature of the PMSM
%
% Syntax
%
% Pfe = ironlosses_PMSM(design, simoptions, vT)
%
% Inputs
%
%   design, simoptions, vT
%
%
% Formula for losses taken from:
% 
% H. Polinder, M. E. C. Damen and F. Gardner, 'Design, modelling and test
% results of the AWS PM linear generator',  European Transactions on
% Electrical Power, 2005, vol 15, pp 245–256, Published online in Wiley
% InterScience (www.interscience.wiley.com). DOI: 10.1002/etep.56
%

    fe = velocity2electricalfreq(vT, design.PoleWidth);
    
    mt = design.Wt * design.ht * design.ls * design.phases ...
         * design.PowerPoles * design.sides * simoptions.evaloptions.ArmatureIronDensity;
    
    my = design.hba * design.Wp * design.ls * design.PowerPoles ...
         * design.sides * simoptions.evaloptions.ArmatureIronDensity;
    
    Pfe = 2 .* design.Pfe0 ...
          .* (mt .* ((design.Ws + design.Wt)/design.Wt).^2 ...
              + my .* (design.Wp ./ (pi * design.hba)).^2) ...
          .* (fe ./ design.fe0) .* (design.BgPeak / design.Bg0).^2;
   
end

function fe = velocity2electricalfreq(v, polewidth)
% converts a linea machine velocity to its electrical frequency at this
% speed
%
% Syntax
%
% fe = velocity2electricalfreq(v, polewidth)
%
% Input
%
%   v - the relative velocity of the two machine parts

    % s = d/t
    % st = d
    % t = d / s
    % f = s / d
    fe = abs(v) ./ (2 * polewidth);

end