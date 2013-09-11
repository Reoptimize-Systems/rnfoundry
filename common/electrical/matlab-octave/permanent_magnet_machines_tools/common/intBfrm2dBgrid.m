function intB = intBfrm2dBgrid(x, y, w, h, coords, X, Y, B, edgemode, scalefac)
% calculates the integral of the flux density over a rectangular
% cross-section from a half-period 2D flux density polynomial at several
% positions
%
% Syntax
% 
% intB = intBfrm2dBgrid(x, y, w, h, coords, X, Y, B)
% intB = intBfrm2dBgrid(x, y, w, h, coords, X, Y, B, edgemode)
% intB = intBfrm2dBgrid(x, y, w, h, coords, X, Y, B, edgemode, scalefac)
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
%  X, Y - (n x p) matrices containing x and y position data in the same
%    format as produced by  meshgrid(xnodes,ynodes)
%
%  B - (n x p x 2) matrix of sets of values of the flux density in the x
%    and y direction at the locations specified in X, Y. B(:,:,1) contains
%    the x-directed flux density at the positions in X,Y and B(:,:,2)
%    contains the y-directed flux density at the positions in X,Y.
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
%    integrals for each of the n coordinates provided in coords.
%
    
    if nargin < 10
        scalefac = 1;
    end
    
    if nargin < 9
        edgemode = true;
    end
    
    intB = zeros(size(coords,1), 2);
    
    if any(diff(coords(:,1)))
        % the x coordinates change, and therefore the region of integration
        % in the x direction changes at each coordinate, we must create a
        % full solution for every coordinate
        
        for i = 1:size(coords,1)

            intB(i,1) = integratehalfperiod2ddata(X, Y, B(:,:,1), ...
                                         x + coords(i,1), ...
                                         y + coords(i,2), ...
                                         x + w + coords(i,1), ...
                                         y + h + coords(i,2), ...
                                         edgemode);


            intB(i,2) = integratehalfperiod2ddata(X, Y, B(:,:,2), ...
                                         x + coords(i,1), ...
                                         y + coords(i,2), ...
                                         x + w + coords(i,1), ...
                                         y + h + coords(i,2), ...
                                         ~edgemode);

        end
    
    else
        % if the x positions are the same we can use the precomputed x line
        % integral slm made in the first call to integratehalfperiod2ddata
        % in each subsequent call as the region of integration in the x
        % direction is not changing
        
        [intB(1,1), intBslm1] = integratehalfperiod2ddata(X, Y, B(:,:,1), ...
                                     x + coords(1,1), ...
                                     y + coords(1,2), ...
                                     x + w + coords(1,1), ...
                                     y + h + coords(1,2), ...
                                     edgemode);


        [intB(1,2), intBslm2] = integratehalfperiod2ddata(X, Y, B(:,:,2), ...
                                     x + coords(1,1), ...
                                     y + coords(1,2), ...
                                     x + w + coords(1,1), ...
                                     y + h + coords(1,2), ...
                                     ~edgemode);
        
        for i = 2:size(coords,1)
            
            intB(i,1) = integratehalfperiod2ddata(intBslm1, ...
                                                 y + coords(i,2), ...
                                                 y + h + coords(i,2));


            intB(i,2) = integratehalfperiod2ddata(intBslm2, ...
                                                 y + coords(i,2), ...
                                                 y + h + coords(i,2));

        end
        
        
    end
    
    % multiply by a scale factor (usually used to denormalise the result)
    intB = intB .* scalefac;
    
end