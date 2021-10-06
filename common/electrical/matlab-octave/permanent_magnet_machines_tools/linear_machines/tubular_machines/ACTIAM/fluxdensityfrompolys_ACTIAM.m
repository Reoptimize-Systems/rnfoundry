function FluxDensity = fluxdensityfrompolys_ACTIAM(WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, Rm, rzCoords, Polynomials)
% Returns the normalised flux density above the translator. This is
% calculated by evaluating the polynomial approximations of the flux
% density at the appropriate positions.
%
% arguments: (input)
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
%   rzCoords - (n x 2) matrix of r and z coordinates at which the
%              flux density is to be evaluated. The r coordinates must be
%              normalised such that r is stated relative to Rm, i.e. an r
%              position 1 cm above a translator of Rm = 10 cm would be
%              restated as r/Rm = 1/10 = 0.1. The z coordinates must be
%              restated as ratios of a pole width in a similar manner where
%              the actual distance from the centre of a steel piece is
%              divided by the total pole width. Values up to z = 1 are
%              allowed.
%
%   polynomials - polynomials to be used to evaluate the flux density, if
%                 not supplied they will be loaded from disc
%
% Arguments: (output)
%
%   FluxDensity - (n x 2) Matrix of Flux Density values at the coordinates
%     supplied in the input. The first column is the axial flux at the
%     given points, the second the radial flux.
%

% Copyright 2007 Richard Crozier and The Institute For Energy Systems at
% The University of Edinburgh


    % First the polynomials which will be used to evaluate the flux
    % density must be loaded into memory if they have not been passed 
    % to the function by it's parent function. These are stored in a
    % (i x 4 x 204) array. The row size varies for each polynomial.
    % The first three columns in each of the 204 polynomials
    % corresponds to a variable describing the machine and the
    % location above the translator surface i.e. Wm/Wp, Wp/Rm and r/Rm
    % where the row values are the power to which variable must
    % be raised. The forth column contains the coefficients by which
    % each term must be multiplied. The first 51 polynomials are close
    % axial polynomials, 52 to 102 the close radial positions, 103 to
    % 153 the far axial polynomials and 154 to 204 the far radial
    % polynomials.
    
    if nargin < 8
        load('BPolynomials_IA.mat');
    end
    
    % Next we should calculate which polynomials must be used at
    % each r,z Coordinate position. GetPolys returns an (n x 8)
    % matrix where columns 1, 2, 5 and 6 are the 3rd dimension
    % indices of p that should be used to evaluate the flux density
    % in the axial and radial directions respectively
    % and the other columns are the corresponding equivalent z/Wp
    % values of the polynomials in order to facilitate interpolation.

    polys = selectpolys_TM(rzCoords);

    % Next the polynomials can be evaluated for the given positions
    % and machine properties. evaluatepolys_ACTIAM returns an (n x 2)
    % matrix, the first column corresponding to the first column of
    % polys, the second, the second column of polys.

    tempFlux = evaluatepolys_ACTIAM(WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, Rm, rzCoords(:,1), polys(:,1:2), Polynomials);
    tempFlux(:,3:4) = evaluatepolys_ACTIAM(WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, Rm, rzCoords(:,1), polys(:,3:4), Polynomials);

    % Next we use the returned values and their relative positions to
    % interpolate the flux density at the actual coordinates.

    FluxDensity = zeros(size(tempFlux,1), 2); % Preallocate for speed

    for n = 1:size(tempFlux,1)

        % First find axial Flux density
        FluxDensity(n,1) = interp1(polys(n,5:6),tempFlux(n,1:2),abs(rzCoords(n,2)),'spline');

        % If the sign of the z coordinate is negative, we must reverse the
        % polarity of the flux
        if rzCoords(n,2) < 0

            FluxDensity(n,1) = FluxDensity(n,1) * -1;

        end

        % Next find radial Flux density
        FluxDensity(n,2) = interp1(polys(n,5:6),tempFlux(n,3:4),abs(rzCoords(n,2)),'spline');

        % If the sign of the z coordinate is negative, we must reverse the
        % polarity of the flux
        if rzCoords(n,2) < 0

            FluxDensity(n,2) = FluxDensity(n,2) * -1;

        end

    end

end