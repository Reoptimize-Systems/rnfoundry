function varargout = AxialLoadCaCoeffs(P, E, I, a, l)
% AxialLoadFCoeffs: calculates Ca1 to Ca6 for the calculation of
% deflections etc. for beams under axial and transverse loading
%
% Input:
%
%   P - axial load at A
%
%   E - young's modulus of beam material
%
%   I - second moment of area of the beam
%
%   a - position of load on beam
%
%   l - length of beam 
%
% Output:
%
%   Can - (1 x 6) matrix of the Ca coefficients for the beam
%
%   k - optional, k value

    k = sqrt(P./(E.*I));
    
    Can = zeros(6, size(k,2));
    
    Can(1,:) = cos(k.*(l - a));
    
    Can(2,:) = sin(k.*(l - a));
    
    Can(3,:) = 1 - cos(k.*(l - a));
    
    Can(4,:) = k.*(l - a) - sin(k.*(l - a));
    
    Can(5,:) = k.^2.*((l - a).^2)./2 - Can(3);
    
    Can(6,:) = k.^3.*((l - a).^3)./6 - Can(4);
    
    varargout{1} = Can;
    varargout{2} = k;
    
end