function J = SolidCircle (JVars)
% caclulates the polar moment of inertail of a solid circular section
% according to roark's formulas for stress and strain
%
% Syntax
%
% J = roark.PolarMoment.SolidCircle (JVars)
%
% Input
%
%  JVars - (n x 1) column matrix of values of R the radii of the solid
%   circles for which the polar moment of inertia is to be calculated
%
% Output
%
%  J - (n x 1) values of polar moment of inertia for each circle
%
%

    J = 0.5 .* pi .* realpow (JVars(:,1), 4);

end