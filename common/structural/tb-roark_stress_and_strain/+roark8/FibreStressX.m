function sigma = FibreStressX (Ivars, Yvars, x, IMethod, beamMethod, E)
% calculates the maximum fibre in a beam about the x axis.
%
% Syntax
%
% sigma = FibreStressX (Ivars, Yvars, x, IMethod, beamMethod)
% sigma = FibreStressX (Ivars, Yvars, x, IMethod, beamMethod, E)
%
% Input
%
%   Ivars - (n x p) matrix of values necessary for caluculating the
%     second moment of inertia according to the method described in
%     'IMethod'. See the appropriate function for details of the
%     required variables, and exact format of Ivars.
%
%   Yvars - (n x p) matrix of values necessary for caluculating the
%     beam moments according to the method described in 'beamMethod'. See
%     the appropriate function for details of the required variables, and
%     exact format of Yvars.
%
%   x - row vector of position values at which the deflection is to be
%     calculated
%
%   beamMethod - string describing the method by which the beam deflection
%     is to be calculated. These should correspond to the appropriate table
%     in Roark's Formulas for Stress & Strain.
%
%   IMethod - string describing the method by which the second moment of
%     inertia is to be calculated. These should correspond to the
%     appropriate table in Roark's Formulas for Stress & Strain. Only
%     required for some cases.
%
%   E - Young's modulus of the beam material, only required in some cases
%
% Output
%
%   sigma - (n x 1) column vector of values of the maximum fibre stress in
%     a beam, at the positions specified in 'x', calulated according to
%     'IMethod' and 'beamMethod'.
%

    if nargin < 6
        Mom = roark8.BeamMomentSuperX (Yvars, x, beamMethod);
    else
        Mom = roark8.BeamMomentSuperX (Yvars, x, beamMethod, Ivars, IMethod, E);
    end
    
    I = roark8.MomentOfInertiaX (IVars, IMethod);
    
    y_c = roark8.DistanceToExtremetiesX (IVars, IMethod);
    
    sigma = -Mom .* y_c ./ I;
    
end