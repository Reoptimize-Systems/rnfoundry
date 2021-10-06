function varargout = AxialLoadFaCoeffs(P, E, I, a, x)
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
%   a - position of load on beam
%
%   x - position along beam 
%
% Output:
%
%   Fan - (6 x n) matrix of the Fa coefficients for each x position
%         specified in x

    k = sqrt(P./(E.*I));
    
    Fan = zeros(6, size(x,2));
    
    Fan(1, x >= a) = cos(k.*(x(x>=a) - a));
    
    Fan(2, x >= a) = sin(k.*(x(x>=a) - a));
    
    Fan(3, x >= a) = (1 - cos(k.*(x(x>=a) - a)));
    
    Fan(4, x >= a) = k.*(x(x>=a) - a) - sin(k.*(x(x>=a) - a));
    
    Fan(5, x >= a) = ((k.^2) .* ((x(x>=a) - a).^2)) ./ 2 - Fan(3,x>=a);
    
    Fan(5, x < a) = 0 - Fan(3,x<a);
    
    Fan(6, x >= a) = ((k.^3) .* ((x(x>=a) - a).^3)) ./ 6 - Fan(4,x>=a);
          
    Fan(6, x < a) = 0 - Fan(4,x<a);
    
    varargout{1} = Fan;
    varargout{2} = k;
    
end