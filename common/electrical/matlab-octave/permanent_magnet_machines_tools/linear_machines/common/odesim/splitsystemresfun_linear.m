function [results, design] = splitsystemresfun_linear(T, Y, design, simoptions, results)
% splitsystemresfun_linear: finalizes the results for a linear machine
% simulation as accumulated by splitodesystemres_linear, using odesplit
% to reduce the memory demands
%
% Syntax
%
% [results, design] = splitsystemresfun_linear(T, Y, design, simoptions, results)
%
%

        % copy over the minimum length of the longer machine part
        design.minLongMemberLength = results.minLongerPartLength;
        % calculate the number of pole widths required to achieve this
        % length
        design.minLongMemberPoles = ceil(design.minLongMemberLength ./ design.PoleWidth);
        % recalculate the length to reflect the quantized nature of the
        % poles
        design.minLongMemberLength = design.minLongMemberPoles * design.PoleWidth;
        % copy over the max experienced relative speed
        design.vRmax = results.vRmax;
        % now do common split system results function tasks
        [results, design] = splitsystemresfun_AM(T, Y, design, simoptions, results);
        
end
