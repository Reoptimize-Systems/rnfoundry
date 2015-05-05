function strandDc = stranddiameter (Dc, Nstrands)
% calculates the required strans diameter to make the same area of wire as
% an unstranded conductor
%
% Syntax
%
% strandDc = stranddiameter (Dc, Nstrands)
%
% Input
%
%   Dc - solid conductor diameter
%
%   Nstrands - desired number of parallel strands to make same area
%
%

    warea = circlearea (Dc./2);
    
    strandarea = warea ./ Nstrands;
    
    strandDc = 2 .* area2radius (strandarea);

end