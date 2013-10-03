function varargout = AxialLoadCCoeffs(P, E, I, l)
% AxialLoadFCoeffs: calculates C1, C2, C3, and C4 for the calculation of
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
%   l -length of beam 
%
% Output:
%
%   Cn - C coefficients C1, C2, C3 and C4 in a (1 x 4) vector
%
%   k - optional, k value

    k = sqrt(P./(E.*I));
    
    Cn = zeros(4,size(k,2));
    
    Cn(1,:) = cos(k.*l);
    Cn(2,:) = sin(k.*l);
    Cn(3,:) = 1 - cos(k.*l);
    Cn(4,:) = (k.*l) - sin(k.*l);
    
    varargout{1} = Cn;
    varargout{2} = k;
end