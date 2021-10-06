function CVV = CVVcoeff(l2, E2, I2, l3, E3, I3)
% CVVcoeff: calculates the CVV coefficient for a frame with it's right side
% fixed at the base for the cases 5 to 12 in Table 4 of Roark's Formulas
% for Stress and Strain
%
% Input:
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
    CVV = (l2.*(l3.^2))./(E2.*I2) + (l3.^3)./(3.*E3.*I3); 

end