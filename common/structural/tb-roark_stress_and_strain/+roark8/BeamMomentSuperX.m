function Mom = BeamMomentSuperX (Yvars, x, beamMethod, Ivars, IMethod, E)
% calculates the moments in a beam using the superposition of two or more
% loads
%
% Syntax
%
% Mom = BeamMomentSuperX (Yvars, x, beamMethod)
% Mom = BeamMomentSuperX (Yvars, x, beamMethod, Ivars, IMethod, E)
%
% Input
%
%   Yvars - If using identical method for all superposition cases:
%
%     An (n x p) matrix of values necessary for caluculating the second
%     moment of inertia according to the method described in 'IMethod'. See
%     the appropriate function for details of the required variables, and
%     exact format of Ivars.
%
%     If using different methods for all superposition cases, i.e. passing
%     a cell aray of strings to beamMethod:
%
%     An (n x p) cell array of row vectors, each containing a (1 x p)
%     vector of the appropriate of values necessary for calculating the
%     beam moments according to the corresponding string in the cell array
%     'beamMethod'. See each approriate function for details of the
%     required variables, and exact format of each vector.
%
%   x - row vector of position values at which the moment is to be
%     calculated
%
%   beamMethod - string or (n x 1) cell array of strings describing the
%     method by which the beam moment is to be calculated. If a single
%     string, the same method will be used for all cases, if a cell array
%     of strings, each string will define the method to be used for
%     evaluating the parameters in Yvars. These will correspond to the
%     appropriate table in Roark's Formulas for Stress & Strain.
%
%   Ivars - An (1 x p) row vector of values necessary for calculating the
%     second moment of inertia according to the method described in
%     'IMethod'. See the appropriate function for details of the required
%     variables, and exact format of Ivars.
%
%   IMethod - string describing the method by which the second moment of
%     inertia is to be calculated. These should correspond to the
%     appropriate table in Roark's Formulas for Stress & Strain.
%
%   E - Young's modulus of the beam material
%
%
% Output
%
%   Def - (n x i) column vector of values of the moment in a beam, at
%     each of the 'i' positions specified in 'x', calculated according to
%     'IMethod' and 'beamMethod'.
%
    if iscellstr(beamMethod) && size(Yvars, 1) == size(beamMethod, 1)

        Mom = zeros (size (x));
        for i = 1:length(beamMethod)
            Mom = Mom + roark8.BeamMomentSuperX (Yvars{i}, x, beamMethod{i}, Ivars, IMethod, E);
        end

    elseif ischar(beamMethod)

        Mom = sum (roark8.BeamMomentX (Yvars, x, beamMethod, Ivars, IMethod, E),1);

    else
        error('beamMethod must be either a single string or cell array of strings of the same size as Yvars')
    end

end

