function J = HollowCircle (JVars)
% caclulates the polar moment of inertail of a hollow circular section
% according to roark's formulas for stress and strain
%
% Syntax
%
% J = roark.PolarMoment.HollowCircle (JVars)
%
% Input
%
%  JVars - (n x 2) matrix of values of ro and ri, the outer and inner radii
%   of the hollow circles for which the polar moment of inertia is to be
%   calculated respectively.
%
% Output
%
%  J - (n x 1) values of polar moment of inertia for each hollow circle
%
%

    J = 0.5 .* pi .* ( realpow (JVars(:,1), 4) - realpow (JVars(:,2), 4) );

end