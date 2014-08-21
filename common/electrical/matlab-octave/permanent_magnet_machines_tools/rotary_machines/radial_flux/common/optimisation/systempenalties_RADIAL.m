function [score, design, simoptions] = systempenalties_RADIAL(design, simoptions, score)
% adds penalties related to the operation of a radial flux generator system
%
% Syntax
%
% [score, design, simoptions] = systempenalties_RADIAL(design, simoptions, score)
%
% 

    [score, design, simoptions] = systempenalties_ROTARY (design, simoptions, score);

end