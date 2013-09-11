function lambda = fluxlinkagefrm2dApoly(x1, y1, x2, y2, w, h, coords, Apoly, nturns, depth, scalefac)
% calculates the flux linkage in a two part coil circuit from a 2d vector
% potential polynomial at several positions
%
% Syntax
% 
% lambda = fluxlinkagefrm2dA(x1, y1, x2, y2, w, h, coords, Apoly, nturns, depth, scalefac)
% 
% Input
%
%  x1 - x position of bottom left corner of first coil cross-section
% 
%  y1 - y position of bottom left corner of first coil cross-section 
% 
%  x2 - x position of bottom left corner of second coil cross-section
% 
%  y2 - y position of bottom left corner of second coil cross-section 
% 
%  w - the width of the coil cross-section (in x direction)
% 
%  h - the height of the coil cross-section (in x direction) 
% 
%  coords - (n X 2) matrix containing a set of positions at which the coil
%    flux linkage is to be evaluated. These coordinated will be added to
%    the x and y positions of the coil corners, for both coil parts.
% 
%  Apoly -  a scalar structure or vector of 2 structures defining a
%    polynomial fitted to the vector potential in the 2D region where the
%    coils are located. It is these polynomials must be fitted to 
%    a region representing one half period of a regiogion containing a flux
%    waveform periodic in the y direction, and at an absolute peak at the
%    top and bottom of the fitted region. If one polynomial is supplied, it
%    is used to evaluate both coil flux integrals. If two polynomials are
%    supplied, the first is used to evaluate the first coil part, and the
%    second used to evaluate the second coil part. The polynomials must be
%    of the same form as produced by 'polyfitn'.
% 
%  nturns - the number of turns of wire in the coil
% 
%  depth -  the depth of the 2d problem
% 
%  scalefac - (optional) a scalar factor by which to scale the flux linkage
%    results. This can be useful when the position data to which the 
%    polynomials have been fitted has been normalised.
%
% Output
%
%  lambda - (n X 1) vector of flux linkage values for the coi at each
%    corresponding position in coords.
%
    
    if nargin < 11
        scalefac = 1;
    end
    
    tol = [1e-8, 1e-8];
    
    intA = zeros(size(coords,1), 2);
    
    % a separate polynomial can be supplied for each coil part, set
    % appropriate indexes into the array of polynomial structures in each
    % case
    if numel(Apoly) == 1
        polyinds = [1, 1];
    elseif numel(Apoly) > 1
        polyinds = [1, 2];
    end
    
    for i = 1:size(coords,1)
        
        intA(i,1) = integratehalfperiodypoly(Apoly(polyinds(1)), ...
                                             x1 + coords(i,1), ...
                                             y1 + coords(i,2), ...
                                             x1 + w + coords(i,1), ...
                                             y1 + h + coords(i,2), ...
                                             tol(1));
        
        intA(i,2) = integratehalfperiodypoly(Apoly(polyinds(2)), ...
                                             x2 + coords(i,1), ...
                                             y2 + coords(i,2), ...
                                             x2 + w + coords(i,1), ...
                                             y2 + h + coords(i,2), ...
                                             tol(2));
    
        tol = intA(i,:) ./ 100;
        
        tol(tol< 1e-9) = 1e-9;

    end
    
    % flux linkage in a 2D planar sim is given by:
    %
    % $$ \frac{N l_s}{S} * ( int_S A_+  - int_S  A_- ) $$
    lambda = nturns .* depth .* (intA(:,1) - intA(:,2)) .* scalefac ./ (w * h); 
    
end