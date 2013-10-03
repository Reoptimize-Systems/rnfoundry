function Mom = GenericMomentDistribLoad(MA, RA, wa, wl, x, l, a, E, I)
% GenericMomentDistribLoad: Calculates the moment in a beam with some
% distributed loading using the generic moment formula for distributed
% loads described in Table 2, row header row 2 on page 102 in 'Roark's
% Formulas Stress & Strain 6th edition'. You are required to supply
% information such as the reaction force and moments etc. See Roark for a
% full description of the inputs below.
%
% Input: 
%   
%   MA - applied couple (moment) at point A
%
%   RA - Reaction force at point A
%
%   wa - unit load at 'a'
%
%   wl - unit load at M_B, the end of the beam
%
%   x - matrix of position values at which the deflection is to be calculated 
%
%   l - length of the beam
%
%   a - distance from M_A at which 'wa' is applied 
%
%   E - Young's modulus of the beam material
%
%   I - second moment of inertia of the beam cross-section
%
% Output:
%
%   Mom - values of the moment
%
    if x < a
        stepfun1 = 0;
        stepfun2 = 0;
    else
        stepfun1 = (x-a).^2;
        stepfun2 = (x-a).^3;
    end
    
    Mom = MA + RA .* x - wa .* stepfun1 ./ 2 ...
          - (wl - wa) .* stepfun2 ./ (24 .* E .* I .* (l - a));
      
end

