function SheathForce = sheathforcefrompolys_ACTIAM(WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, Rm, div, Polynomial)
% Evaluates the force on the outer sheath of an ACTIAM from fitted
% polynomials
%

    if nargin < 8
        load('SheathForcePoly_IA.mat');
    end

    % this is the total per-pole force due to the sheath
    SheathForce = evalforcepoly(WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, Rm, Polynomial);
    
    % This is then divided by the desired number of divisions
    SheathForce = SheathForce ./ div;
    
end

function SheathForce = evalforcepoly(WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, Rm, p)
% Evaluates the polynomials for the positions and specifications of a
% particular machine.
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
%   p - the array of polynomials produced from approximations of FEA
%       results and stored in the Polynomials matrix
%
% Arguments: (output)
%        NormSheathForce - (n x k) column vector of values of the
%        normalised sheath Force

    % First we can use the Wm/Wp and Wp/Rm values to reduce the calculation
    % necessary for each position by restating the coefficients of the
    % polnomial by inputting these values as constants.

    vars = repmat([WmVWp WpVRm RoVRm RaVRo RsoVRm Rm],size(p,1),1);
    
    % The results of the above operations are then summed to give 
    % the force
    SheathForce = sum(prod([p(:,1) realpow(vars,p(:,2:end))],2),1);


end
