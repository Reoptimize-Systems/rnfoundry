function Mom = BeamMomentX (Yvars, x, beamMethod, Ivars, IMethod, E)
% calculates the moment in a beam about the x axis.
%
% Syntax
%
% Mom = BeamMomentX (Yvars, x, beamMethod, Ivars, IMethod, E)
%
% Input
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
% The following arguments are only required for certain cases such as beams
% undergoing both lateral and axial laoding.
%
%   Ivars - (n x p) matrix of values necessary for caluculating the
%     second moment of inertia according to the method described in
%     'IMethod'. See the appropriate function for details of the
%     required variables, and exact format of Ivars. Only required for some
%     cases.
%
%   IMethod - string describing the method by which the second moment of
%     inertia is to be calculated. These should correspond to the
%     appropriate table in Roark's Formulas for Stress & Strain. Only
%     required for some cases.
%
%   E - Young's modulus of the beam material, only required in some cases
%
% Output:
%
%   Mom - (n x 1) column vector of values of the moment in a beam, at
%     the positions specified in 'x', calulated according to 'IMethod'
%     and 'beamMethod'.
%

    switch beamMethod

        case {'8.1.1a', 'LURFP', 'lurfp'}
            % Left end free (unsupported), right end fixed,
            % point load
            Mom = roark.Beams.ConcLoad.LURF.BendingMoment (Yvars, x);
            
        case {'8.1.1e', 'LSRSP', 'lsrsp'}
            % R6 T3.1e
            % Left end Simply Supported, right end simply supported,
            % point load
            Mom = roark.Beams.ConcLoad.LSRS.BendingMoment (Yvars, x);
            
        case {'8.1.1f', 'LGRSP', 'lgrsp'}
            % R6 T3.1f
            % Left end guided, right end simply supported,
            % point load
            Mom = roark.Beams.ConcLoad.LGRS.BendingMoment (Yvars, x);
            
        case {'8.8.2e', 'LSRSDAL', 'lsrsdal'}
            % R6 T10.2e
            % Left end Simply Supported, right end simply supported,
            % distributed force and axial load
            I = roark8.MomentOfInertiaX (Ivars, IMethod);
            
            Mom = Table10r2eMom (Yvars, E, I, x);
            
        case {'8.1.3a', 'LURFM', 'lurfm'}
            % Left end free, right end fixed,
            % applied moment
            Mom = roark.Beams.ConcMoment.LURF.BendingMoment (Yvars, x);
            
        otherwise
            
            Mom = feval(beamMethod, Yvars, x);

    end

end