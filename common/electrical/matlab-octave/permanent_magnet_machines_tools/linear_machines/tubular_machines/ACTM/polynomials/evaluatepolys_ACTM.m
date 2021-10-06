function FluxDensity = evaluatepolys_ACTM(WmVWp, WpVRm, RsoVRm, Rm, rCoords, polys, p)
% evaluatepolys_ACTM: a function to evaluate the polynomials for the
% positions and specifications of a particular machine.
%
% Arguments: (input)
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
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
%           they are storedin the multidimensional array 'p'. In effect 
%           they are the location in the third dimension of p where the
%           matrix corresponding to the appropriate polynomial is stored.
%
%   p - the array of polynomials produced from approximations of FEA
%       results and stored in the Polynomials.mat
%       file.
%
% Arguments: (output)
%        FluxDensity - (n x k) column vector of Flux Density values at the
%        'n' r coordinates supplied in rCoords for each of the 'k'
%        polynomials to be evaluated. Typically this will be two
%        polynomials, i.e. k = 2, as fluxdensityfrompolys_ACTM, for which this program
%        had been designed interpolates between two polynomials in order to
%        find the flux  

    vars = repmat([WmVWp WpVRm RsoVRm Rm 0],size(p,1),1);
    
    FluxDensity = zeros(size(rCoords,1), size(polys, 2)); % Preallocate for speed

    for n = 1:size(rCoords,1)
        % For each r Coordinate each polynomial on the corresponding row
        % of polys must be evaluated.

        for k = 1:size(polys, 2)
            % We raise the values of 'r' to the powers in the second column
            % of the reduced polynomials then multiply these values by the
            % coefficients in the first column.
            vars(:,end) = rCoords(n,1);
            % The results of the above operations are then summed to give 
            % the flux density at the point.
            FluxDensity(n,k) = sum(prod([p(:,1,polys(n,k)) realpow(vars,p(:,2:end,polys(n,k)))],2),1);
        end

    end

end
