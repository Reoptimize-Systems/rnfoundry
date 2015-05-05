function J = SolidCircle (JVars, T)
% calculates the maximum tosional shear stress in a section according to
% roark's formulas for stress and strain
%
% Syntax
%
% J = roark.TorsionalShearStress.SolidCircle (JVars)
%
% Input
%
%  JVars - (n x 1) column matrix of values of R the radii of the solid
%   circles for which the stress is to be calculated
%
% Output
%
%  J - (n x 1) values of maximum shear stress for each circle
%
%

    J = 2 .* T ./ ( pi .* realpow (JVars(:,1), 3) );

end