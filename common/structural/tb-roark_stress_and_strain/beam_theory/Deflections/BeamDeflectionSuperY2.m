function Def = BeamDeflectionSuperY2(Yvars, Ivars, E, x, IMethod, beamMethod)
% function: BeamDeflectionSuper
%
% Function for calculating the deflection of a beam using the superposition
% of two or more loads
%
% Input:
%
%   Ivars - An (1 x p) row vector of values necessary for caluculating the
%           second moment of inertia according to the method described in 
%           'IMethod'. See the appropriate function for details of the
%           required variables, and exact format of Ivars.
%
%   Yvars - If using identical method for all superposition cases:
%
%           An (n x p) matrix of values necessary for caluculating the
%           second moment of inertia according to the method described in 
%           'IMethod'. See the appropriate function for details of the
%           required variables, and exact format of Ivars.
%           
%           If using different methods for all superposition cases, i.e. 
%           passing a cell aray of strings to beamMethod:
%
%           An (n x p) cell array of row vectors, each containing a (1 x p) 
%           vector of the appropriate of values necessary for calculating the
%           beam deflection according to the corresponding string in the cell 
%           array 'beamMethod'. See each approriate function for details of the
%           required variables, and exact format of each vector.
%
%   E - Young's modulus of the beam material
%
%   x - row vector of position values at which the deflection is to be
%       calculated 
%
%   IMethod - string describing the method by which the second moment of 
%            inertia is to be calculated. These should correspond to the 
%            appropriate table in Roark's Formulas for Stress & Strain.
%
%   beamMethod - string or (n x 1) cell array of strings describing the 
%                method by which the beam deflection is to be calculated.
%                If a single string, the same method will be used for all
%                cases, if a cell array of strings, each string will define
%                the method to be used for evaluating the parameters in
%                Yvars. These will correspond to the appropriate table in 
%                Roark's Formulas for Stress & Strain. 
% Output:
%
%   Def - (n x i) column vector of values of the deflection of a beam, at
%         each of the 'i' positions specified in 'x', calulated according 
%         to 'IMethod' and 'beamMethod'.
%    
    if iscellstr(beamMethod) && size(Yvars, 1) == size(beamMethod, 1)
        
        Def = 0;
        for i = 1:length(beamMethod)
            Def = Def + BeamDeflectionY2(Ivars, Yvars{i,1}(1,:), E, x, IMethod, char(beamMethod(i,1)));
        end

    elseif ischar(beamMethod)

        Def = sum(BeamDeflectionY2(Ivars, Yvars, E, x, IMethod, beamMethod),1);

    else  
        error('beamMethod must be either a single string or cell array of strings of the same size as Yvars')
    end

end

