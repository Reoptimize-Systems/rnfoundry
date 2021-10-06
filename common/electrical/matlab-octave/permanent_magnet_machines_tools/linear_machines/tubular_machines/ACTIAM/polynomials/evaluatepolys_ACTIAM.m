function FluxDensity = evaluatepolys_ACTIAM(WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, Rm, rCoords, polys, p)
% evaluate the polynomials for the positions and specifications of an
% slotless tubular permanent magnet macine.
%
% Arguments: (input)
%
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RoVRm - scalar value of Ro/Rm Ratio for machine to be evaluated, in
%           order to define the coil height
%
%   RaVRo - scalar value of Ra/Ro Ratio for machine to be evaluated, in
%           order to define the coil height
%
%   RsoVRm - sclar value of Rso/Rm, the ratio of the shaft outer diameter
%            to the translator radius
%
%   Rm - Radius of translator in m
%
%   rCoords - column vector or scalar of r positions at which the
%             evaluation of a polynomial (or polynomials) is to be
%             performed.
%
%   polys - (n x k) column vector or matrix of polynomial numbers. 
%           These are the numbers of the polynomials to be evaluated as 
%           they are stored in the multidimensional array 'p'. In effect 
%           they are the location in the third dimension of p where the
%           matrix corresponding to the appropriate polynomial is stored.
%
%   p - the array of polynomials produced from approximations of FEA
%       results and stored in the Polynomials matrix
%
% Arguments: (output)
%
%   FluxDensity - (n x k) column vector of Flux Density values at the
%                 'n' r coordinates supplied in rCoords for each of the 'k'
%                 polynomials to be evaluated. 
%

    % First we can use the Wm/Wp and Wp/Rm values to reduce the calculation
    % necessary for each position by restating the coefficients of the
    % polnomial by inputting these values as constants.

    vars = repmat([WmVWp WpVRm RoVRm RaVRo RsoVRm Rm 0],size(p,1),1);
    
    FluxDensity = zeros(size(rCoords,1), size(polys, 2)); % Preallocate for speed
    
    %temp = zeros(816,1);

    for n = 1:size(rCoords,1)

        % For each r Coordinate each polynomial on the corresponding row
        % of polys must be evaluated.
        vars(:,end) = rCoords(n,1);
        
        for k = 1:size(polys, 2)

            % We raise the values of 'r' to the powers in the second column
            % of the reduced polynomials then multiply these values by the
            % coefficients in the first colum.
            FluxDensity(n,k) = sum(prod([p(:,1,polys(n,k)) realpow(vars,p(:,2:end,polys(n,k)))],2),1);
            
        end

    end

end
