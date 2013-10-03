function varargout = AxialLoadFCoeffs(P, E, I, x)
% AxialLoadFCoeffs: calculates F1, F2, F3, and F4 for the calculation of
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
%   x - position along beam 
%
% Output:
%
%   

    k = sqrt(P./(E.*I));
    
    Fn = cos(k.*x);
    Fn(2,:) = sin(k.*x);
    Fn(3,:) = 1 - cos(k.*x);
    Fn(4,:) = (k.*x) - sin(k.*x);
    
    varargout{1} = Fn;
    varargout{2} = k;
    
end