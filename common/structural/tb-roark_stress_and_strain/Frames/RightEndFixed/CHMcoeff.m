function CHM = CHMcoeff(l1, E1, I1, l2, E2, I2, l3, E3, I3)
% CHMcoeff: calculates the CHM coefficient for a frame with it's right side
% fixed at the base for the cases 5 to 12 in Table 4 of Roark's Formulas
% for Stress and Strain
%
% Input:
%
%   l1 - the length of the right vertical member
%
%   E1 - Young's modulus of the right vertical member
%
%   I1 - Second moment of area of the right vertical member
%
%   l2 - the length of the left vertical member
%
%   E2 - Young's modulus of the lef vertical member
%
%   I2 - Second moment of area of the left vertical member
%
%   l3 - the length of the horizontal member
%
%   E3 - Young's modulus of the horizontal member
%
%   I3 - Second moment of area of the horizontal member
%
    CHM = (l1.^2) ./ (2 .* E1 .* I1) + ((l2./(2.*E2.*I2)).*(2.*l1-l2)) + (l1.*l3./(E3.*I3)); 

end