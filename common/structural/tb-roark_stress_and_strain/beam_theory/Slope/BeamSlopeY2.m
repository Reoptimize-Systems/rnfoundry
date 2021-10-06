function Slope = BeamSlopeY2(Ivars, Yvars, E, x, IMethod, beamMethod)
% function: BeamSlopeY2
%
% function for calculating the slope of a beam.
%
% Input:
%
%   Ivars - (n x p) matrix of values necessary for caluculating the
%          second moment of inertia according to the method described in
%          'IMethod'. See the appropriate function for details of the
%          required variables, and exact format of Ivars.
%
%   Yvars - (n x p) matrix of values necessary for caluculating the
%          beam deflection according to the method described in
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
%   Slope - (n x 1) column vector of values of the slope of a beam, at
%           the positions specified in 'x', calulated according to
%           'IMethod' and 'beamMethod'.
%

    I = MomentOfInertiaY2(Ivars, IMethod);

    switch beamMethod

        case '3.1f'
            % Left end guided, right end simply supported (Cantilever),
            % point load
            Slope = Table3r1fSlope(Yvars, E, I, x);

        otherwise

            Slope = feval(beamMethod, Yvars, E, I, x);

    end

end