function LFH = Table4r5bLFH(Yvars, E1, I1, E2, I2, E3, I3)
% Input: 
%   
%   Yvars - (n x 1) column vector of values of R, the radius of the
%          circular cross-section:
%          Yvars(:,1) - wa, unit load at M_A, the left end of the frame
%          Yvars(:,2) - wb, unit load at M_B, the right end of the frame
%          Yvars(:,3) - l1, length of left vertical member
%          Yvars(:,4) - l2, length of right vertical member 
%          Yvars(:,5) - l3, length of the horizontal beam
%
%   E1 - Young's modulus of the beam material
%
%   I1 - second moment of inertia of the beam cross-section
%
%   E2 - Young's modulus of the beam material
%
%   I2 - second moment of inertia of the beam cross-section
%
%   E3 - Young's modulus of the beam material
%
%   I3 - second moment of inertia of the beam cross-section

    wa = Yvars(:,1);
    wb = Yvars(:,1);
    l1 = Yvars(:,1);
    l2 = Yvars(:,1);
    l3 = Yvars(:,1);  

    LFH = wa.*(l2.*l3.^2.*(2.*l1 - l2)./(4.*E2.*I2) + l1.*l3.^3./(6.*E3.*I3))...
        + (wb - wa).*(l2.*l3.^2.*(2.*l1 - l2)./(12.*E2.*I2) + l1.*l3.^3./(24.*E3.*I3));
    
end