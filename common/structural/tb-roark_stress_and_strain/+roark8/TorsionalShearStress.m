function Tmax = TorsionalShearStress (vars, JMethod, T)
% calculates the maximum torsional shear stress in a section
%
% Syntax
%
% Tmax = roark8.TorsionalShearStress (vars, JMethod, T)
%
% Input:
%
%   vars - (n x p) column vector of values necessary for calculating the
%     second moment of inertia according to the method described in
%     'method'. See the appropriate function for details of the required
%     variables, and exact format of vars.
%
%   JMethod - string describing the method by which the second moment of 
%     inertia is to be calculated. These should correspond to the
%     appropriate table in Roark's Formulas for Stress & Strain, 8th
%     edition. If not one of the listed functions, MomentOfArea will
%     interpret this as a function name and attmpt to evaluate the named
%     function using IVars as the input parameters. Aliases with human
%     understandable names are also available for most sections. Currently,
%     the following sections are available:
%
%           method string options             IVars
%
%     '10.1.1', 'SolidCircle','Circle'     [ R ]
%     '10.1.10', 'Annulus', 'HollowCircle' [ Ro, Ri ]
%
%   T - (n x p) column vector of values of T the moment applied to the
%     section for each section.
%
% Output:
%
%   Tmax - (n x 1) column vector of values of I, the second moment of area of 
%     the cross-section about the x axis, calulated according to 'IMethod'
%
%
% See also: roark.TorsionalShearStress
%

    switch JMethod

        case {'10.1.1', 'SolidCircle', 'Circle'}
            % Solid Circle
            Tmax = roark.TorsionalShearStress.SolidCircle (vars, T);

        case {'10.1.10', 'Annulus', 'HollowCircle'}
            % Hollow Circle
            Tmax = roark.TorsionalShearStress.HollowCircle (vars, T);

        otherwise

            Tmax = feval(JMethod, vars, T);
            
    end

end