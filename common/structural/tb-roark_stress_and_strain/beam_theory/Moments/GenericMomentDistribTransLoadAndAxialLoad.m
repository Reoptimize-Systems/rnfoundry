function Mom = GenericMomentDistribTransLoadAndAxialLoad(thetaA, MA, RA, P, wa, wl, E, I, l, a, x)
% GenericYMomentDistribTransLoadAndAxialLoad: calculates the moment in a
% beam with distributed transverse and concentrated axial loading
% 
% Calculates the moment of a beam with some distributed loading  and a
% concentrated axial load at end A using the generic moment formula for
% distributed loads described in Table 10, row header row 2 on page 164 in
% 'Roark's Formulas Stress & Strain 6th edition'. You are required to
% supply information such as the reaction force and moments etc. See Roark
% for a full description of the inputs below.
%
% Input: 
%   
%   thetaA - externally created concentrated angular displacement at point
%            A in radians
%
%   MA - applied couple (moment) at point A
%
%   RA - Reaction force at point A
%
%   P - axial load applied at A
%
%   wa - unit load at 'a'
%
%   wl - unit load at M_B, the end of the beam
%
%   E - Young's modulus of the beam material
%
%   I - second moment of inertia of the beam cross-section
%
%   l - length of the beam
%
%   a - distance from M_A at which 'wa' is applied 
%
%   x - vector of position values at which the deflection is to be calculated 
%
% Output:
%
%   Mom - values of the deflection 
%
    [Fn, k] = AxialLoadFCoeffs(P, E, I, x);
    
    Fan = AxialLoadFaCoeffs(P, E, I, a, x);

    Mom = MA.*Fn(1,:) + RA.*Fn(2,:)./k - thetaA .* P .* Fn(2,:) ./ k - wa .* Fan(3,:) ./ k.^2 - (wl - wa) .* Fan(4,:) ./ (k.^3 .* (l - a));
    
end