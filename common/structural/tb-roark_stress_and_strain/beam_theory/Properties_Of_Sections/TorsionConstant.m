function K = TorsionConstant(IVars, IMethod)
% function: TorsionConstant
%
% A function for calculating the torsion constant K for a beam
% cross-section. K is a factor dependent on the form and dimensions of the
% cross-section. It is used in the formulas:
%
%     T = (theta / l) * K G      and     theta = T l / KG
%
% where T is the twisting moment, l is the length of the member, theta is
% the angel of twist in radians and G is modulus of rigidity of the
% material.
%
% Input:
%
%   vars - (n x p) column vector of values necessary for calculating the
%          torsion constant according to the method described in 
%          'method'. See the appropriate function for details of the
%          required variables, and exact format of vars.
%
%   method - string describing the method by which the torsion constant is 
%            to be calculated. These should correspond to the  appropriate
%            table in Roark's Formulas for Stress & Strain. If not one of
%            the listed functions, MomentOfArea will interpret this as a
%            function name and attmpt to evaluate the named function using
%            Ivars as the input parameters.
% 
% Output:
%
%   K - (n x 1) column vector of values of K, the torsion constant
%   of the cross-section, calulated according to 'method'
%
    switch IMethod
        
        case '20.4'
            % Solid Square Section
            K = Table20r3K(IVars);
            
        case '20.4'
            % Solid Rectangular Section
            K = Table20r4K(IVars);
            
        case '20.16'
            % Hollow Rectangle
            K = Table20r16K(IVars);
        
        case '1.6'
            % I-Beam
            I = Table1r6I_1(IVars);

        case '1.14'
            % Solid Circle
            I = Table1r14I_1(IVars);

        case '1.15'
            % Hollow Circle
            I = Table1r15I_1(IVars);

        case '1.16'
            % Hollow Circle with very thin annulus
            I = Table1r16I_1(IVars);

        otherwise

            I = feval(IMethod, IVars);
            
    end

end