function J = HollowCircle (JVars)
% calculates the maximum tosional shear stress in a section according to
% roark's formulas for stress and strain
%
% Syntax
%
% J = roark.TorsionalShearStress.HollowCircle (JVars)
%
% Input
%
%  JVars - (n x 1) column matrix of values of R the radii of the hollow
%   circles for which the stress is to be calculated
%
% Output
%
%  J - (n x 1) values of maximum shear stress for each hollow circle
%
%

    J = (2 .* T .* JVars(:,1)) ...
        ./ (pi .* ( realpow (JVars(:,1), 4) - realpow (JVars(:,2), 4) ));

end