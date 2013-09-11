function intB = intBfrm2dBpoly(x, y, w, h, coords, Bpolys, edgemode, scalefac)
% calculates the integral of the flux density over a rectangular
% cross-section from a half-period 2D flux density polynomial at several
% positions
%
% Syntax
% 
% intB = intBfrm2dBpoly(x, y, w, h, coords, Bpolys)
% intB = intBfrm2dBpoly(x, y, w, h, coords, Bpolys, edgemode)
% intB = intBfrm2dBpoly(x, y, w, h, coords, Bpolys, edgemode, scalefac)
% 
% Input
%
%  x - x position of bottom left corner of coil cross-section
% 
%  y - y position of bottom left corner of coil cross-section 
% 
%  w - the width of the coil cross-section (in x direction)
% 
%  h - the height of the coil cross-section (in x direction) 
% 
%  coords - (n X 2) matrix containing a set of positions at which the coil
%    flux linkage is to be evaluated. These coordinated will be added to
%    the x and y positions of the coil corners, for both coil parts.
% 
%  Bpolys -  a vector of 2 structures defining polynomials fitted to the x
%    and y directed flux density in the 2D region where the integral is to
%    be performed. These polynomials must be fitted to a region
%    representing one half period of a region containing a flux waveform
%    periodic in the y direction. The first polynomial is used to evaluate
%    the x directed flux integral, and the second used to evaluate the
%    y-directed flux integral. The polynomials must be of the same form as
%    produced by 'polyfitn'.
% 
%  edgemode - (optional) provides information as to the type of function to be
%    integrated, if edgemode evaluates to true, it is assumed the
%    first polynomial is at a maximum at the edges of the domain in the y
%    direction, if false, it is assumed the first polynomial to be zero-crossing at the
%    edges. Default is true for a maximum at the edges.
%
%  scalefac - (optional) a scalar factor by which to scale the flux density
%    integral results. This can be useful when the position data to which
%    the polynomials have been fitted has been normalised.
%
% Output
%
%  intB - (n X 2) matrix containing the x and y directed flux density
%    integrals for each of n coordinates.
%
    
    if nargin < 8
        scalefac = 1;
    end
    
    if nargin < 7
        edgemode = true;
    end
    
    tol = [1e-8, 1e-8];
    
    intB = zeros(size(coords,1), 2);
    
    for i = 1:size(coords,1)
        
        intB(i,1) = integratehalfperiodypoly(Bpolys(1), ...
                                             x + coords(i,1), ...
                                             y + coords(i,2), ...
                                             x + w + coords(i,1), ...
                                             y + h + coords(i,2), ...
                                             tol(1), ...
                                             edgemode);
        
        intB(i,2) = integratehalfperiodypoly(Bpolys(2), ...
                                             x + coords(i,1), ...
                                             y + coords(i,2), ...
                                             x + w + coords(i,1), ...
                                             y + h + coords(i,2), ...
                                             tol(2), ...
                                             ~edgemode);
    
        tol = intB(i,:) ./ 100;
        
        tol(tol< 1e-9) = 1e-9;

    end
    
    % multiply by a scale factor (usually used to denormalise the result)
    intB = intB .* scalefac;
    
end