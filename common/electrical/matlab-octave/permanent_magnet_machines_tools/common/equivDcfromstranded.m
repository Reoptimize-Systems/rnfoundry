function Dc = equivDcfromstranded (strandDc, Nstrands)
% calculates the required wire diameter to make the same area of wire as
% an stranded conductor
%
% Syntax
%
% Dc = equivDcfromstranded (strandDc, Nstrands)
%
% Input
%
%   strandDc - conductor strand diameter
%
%   Nstrands - desired number of parallel strands to make same area
%
%

    warea = Nstrands .* circlearea (strandDc./2);
    
    Dc = 2 .* area2radius (warea);

end