function Mom = BeamMomentY2(Ivars, Yvars, E, x, IMethod, beamMethod)
% function: BeamMomentY1
%
% function for calculating the moment in a beam about the y1 axis.
%
% Input:
%
%   Ivars - (n x p) matrix of values necessary for caluculating the
%          second moment of inertia according to the method described in
%          'IMethod'. See the appropriate function for details of the
%          required variables, and exact format of Ivars.
%
%   Yvars - (n x p) matrix of values necessary for caluculating the
%          beam moments according to the method described in
%          'beamMethod'. See the appropriate function for details of the
%          required variables, and exact format of Yvars.
%
%   E - Young's modulus of the beam material
%
%   x - row vector of position values at which the deflection is to be
%   calculated
%
%   IMethod - string describing the method by which the second moment of
%            inertia is to be calculated. These should correspond to the
%            appropriate table in Roark's Formulas for Stress & Strain.
%
%   beamMethod - string describing the method by which the beam deflection
%                is to be calculated. These should correspond to the
%                appropriate table in Roark's Formulas for Stress & Strain.
%
% Output:
%
%   Def - (n x 1) column vector of values of the deflection of a beam, at
%         the positions specified in 'x', calulated according to 'IMethod'
%          and 'beamMethod'.
%

    I = MomentOfInertiaY2(Ivars, IMethod);

    switch beamMethod

        case '3.1f'
            % Left end guided, right end simply supported (Cantilever),
            % point load
            Mom = Table3r1fMom(Yvars, E, I, x);
%
%         case '3.2d'
%             % Left end fixed, right end fixed, distributed force
%             Mom = Table3r2dMom(Yvars, E, I, x);
%
%         case '3.2e'
%             % Left end Simply Supported, right end simply supported,
%             % distributed force
%             Mom = Table3r2eMom(Yvars, E, I, x);
%
%         case '3.4d'
%             % Left end fixed, right end fixed, externally created angular
%             % displacement
%             Mom = Table3r4dMom(Yvars, E, I, x);

        case '10.2e'
            % Left end Simply Supported, right end simply supported,
            % distributed force and axial load
            Mom = Table10r2eMom(Yvars, E, I, x);

%         case '10.2f'
%             % Left end guided, right end simply supported,
%             % distributed force and axial load
%             Mom = Table10r2fMom(Yvars, E, I, x);

        otherwise

            Mom = feval(beamMethod, Yvars, E, I, x);

    end

end