function RA = Table3r2dRA(wa, wl, l ,a)
% Table3r2dMA: Calculates the reaction force at end A for a beam with its
% left end fixed and its right end fixed, as calculated in 'Roark's
% Formulas Stress & Strain 6th edition' in table 3, page 104 row 2d.
%
% Input: 
%
%   wa - unit load at 'a'
%
%   wl - unit load at M_B, the end of the beam
%
%   l - length of the beam
%
%   a - distance from M_A at which 'wa' is applied 
%
% Output:
%
%   Mom - values of MA
%

    % Calculate the reaction force at A 
    RA = (wa ./ ((l.^3).*2)) .* ((l-a).^3) .* (l+a) ...
         + (((wl-wa) ./ (20 .* (l.^3)))  .* ((l-a).^3) .* ((3.*l) + (2.*a)) );
    
end