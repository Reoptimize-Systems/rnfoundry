function Def = BeamDeflectionSuperY (Yvars, Ivars, E, x, IMethod, beamMethod)
% Calculates the deflection of a beam in the Y axis using the superposition
% of two or more load cases
%
% Syntax
%
% Def = BeamDeflectionSuperY (Yvars, Ivars, E, x, IMethod, beamMethod)
%
% Input:
%
%   Ivars - An (1 x p) row vector of values necessary for caluculating the
%     second moment of inertia according to the method described in
%     'IMethod'. See the appropriate function for details of the required
%     variables, and exact format of Ivars.
%
%   Yvars - If using identical beamMethod for all superposition cases:
%
%     An (n x p) matrix of values necessary for calculating the beam
%     deflection according to the method described in 'beamMethod'. See the
%     appropriate function for details of the required variables, and exact
%     format of Yvars.
%           
%     If using different methods for all superposition cases, i.e. passing
%     a cell aray of strings to beamMethod:
%
%     An (n x p) cell array of row vectors, each containing a (1 x p)
%     vector of the appropriate of values necessary for calculating the
%     beam deflection according to the corresponding string in the cell
%     array 'beamMethod'. See each approriate function for details of the
%     required variables, and exact format of each vector.
%
%   E - Young's modulus of the beam material
%
%   x - row vector of position values at which the deflection is to be
%     calculated 
%
%   IMethod - string describing the method by which the second moment of 
%     inertia is to be calculated. These should correspond to the
%     appropriate table in Roark's Formulas for Stress & Strain. If not one
%     of the listed functions, MomentOfArea will interpret this as a
%     function name and attmpt to evaluate the named function using IVars
%     as the input parameters. Aliases with human understandable names are
%     also available for most sections. Currently, the following sections
%     are available:
%
%          method string options            IVars
%
%     'A1.2', 'Rectangular'           [ b, d ]
%     'A1.3', 'BoxSection'            [ b, d, di, bi ]
%     'A1.6', 'IBeam'                 [ b, t, tw, d ]
%     'A1.15', 'SolidCircle'          [ R ]
%     'A1.16', 'Annulus'              [ R, Ri ]
%     'A1.17', 'ThinAnnulus'          [ R, t ]
%     'A1.21', 'HollowCircleSector',
%           'AnnularSector'           [ R, t, alpha ]
%     'polygon'                        n x 2 x p matrix 
%
%   beamMethod - string or (n x 1) cell array of strings describing the 
%     method by which the beam deflection is to be calculated. If a single
%     string, the same method will be used for all cases, if a cell array
%     of strings, each string will define the method to be used for
%     evaluating the parameters in Yvars. These will correspond to the
%     appropriate table in Roark's Formulas for Stress & Strain, 8th
%     edition.
%
% Output:
%
%   Def - (n x i) column vector of values of the deflection of a beam, at
%     each of the 'i' positions specified in 'x', calulated according to
%     'IMethod' and 'beamMethod'.
%    

    if iscellstr(beamMethod) && size(Yvars, 1) == numel(beamMethod)
        
        Def = zeros (size (x));
        for i = 1:length(beamMethod)
            Def = Def + roark8.BeamDeflectionSuperY (Ivars, Yvars{i}, E, x, IMethod, beamMethod{i});
        end

    elseif ischar(beamMethod)

        Def = sum (roark8.BeamDeflectionY (Ivars, Yvars, E, x, IMethod, beamMethod),1);

    else  
        error('beamMethod must be either a single string or cell array of strings of the same size as Yvars')
    end

end

