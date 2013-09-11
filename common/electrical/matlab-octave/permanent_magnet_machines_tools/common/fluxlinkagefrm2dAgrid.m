function lambda = fluxlinkagefrm2dAgrid(x1, y1, x2, y2, w, h, coords, X, Y, A, nturns, depth)
% calculates the flux linkage in a two part coil circuit from a 2d vector
% potential polynomial at several positions
%
% Syntax
% 
% lambda = fluxlinkagefrm2dAgrid(x1, y1, x2, y2, w, h, coords, X, Y, A, nturns, depth, scalefac)
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
%  X, Y - (n x p) matrices containing x and y position data in the same
%    format as produced by  meshgrid(xnodes,ynodes)
%
%  A - (n x p x k) matrix of k sets of values of the vetor potential at the
%  locations specified in X, Y, k = 1 or 2. If k == 1 the same data is used
%  to integrate for the two coil regions. If k == 2, A(:,:,1) is used for
%  the first coil part and A(:,:,2) is used for the second part.
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
    
    intA = zeros(size(coords,1), 2);
    
    % a separate data grid can be supplied for each coil part, set
    % appropriate indexes into the array of data matrices in each
    % case
    if size(A,3) == 1
        sliceinds = [1, 1];
    elseif size(A,3) > 1
        sliceinds = [1, 2];
    end
    
    if any(diff(coords(:,1)))
        % the x coordinates change, and therefore the region of integration
        % in the x direction changes at each coordinate, we must create a
        % full solution for every coordinate
        
        for i = 1:size(coords,1)

            intA(i,1) = integratehalfperiod2ddata(X, Y, A(:,:,sliceinds(1)), ...
                                         x1 + coords(i,1), ...
                                         y1 + coords(i,2), ...
                                         x1 + w + coords(i,1), ...
                                         y1 + h + coords(i,2), ...
                                         'max');


            intA(i,2) = integratehalfperiod2ddata(X, Y, A(:,:,sliceinds(2)), ...
                                         x2 + coords(i,1), ...
                                         y2 + coords(i,2), ...
                                         x2 + w + coords(i,1), ...
                                         y2 + h + coords(i,2), ...
                                         'max');

        end
    
    else
        % if the x positions are the same we can use the precomputed x line
        % integral slm made in the first call to integratehalfperiod2ddata
        % in each subsequent call as the region of integration in the x
        % direction is not changing
        
        [intA(1,1), intAslm1] = integratehalfperiod2ddata(X, Y, A(:,:,sliceinds(1)), ...
                                     x1 + coords(1,1), ...
                                     y1 + coords(1,2), ...
                                     x1 + w + coords(1,1), ...
                                     y1 + h + coords(1,2), ...
                                     'max');


        [intA(1,2), intAslm2] = integratehalfperiod2ddata(X, Y, A(:,:,sliceinds(2)), ...
                                     x2 + coords(1,1), ...
                                     y2 + coords(1,2), ...
                                     x2 + w + coords(1,1), ...
                                     y2 + h + coords(1,2), ...
                                     'max');
        
        for i = 2:size(coords,1)
            
            intA(i,1) = integratehalfperiod2ddata(intAslm1, ...
                                                 y1 + coords(i,2), ...
                                                 y1 + h + coords(i,2));


            intA(i,2) = integratehalfperiod2ddata(intAslm2, ...
                                                 y2 + coords(i,2), ...
                                                 y2 + h + coords(i,2));

        end
        
        
    end
    
    
    % flux linkage in a 2D planar sim is given by:
    %
    % $$ \frac{N l_s}{S} * ( int_S A_+  - int_S  A_- ) $$
    lambda = nturns .* depth .* (intA(:,1) - intA(:,2)) ./ (w * h); 
    
end