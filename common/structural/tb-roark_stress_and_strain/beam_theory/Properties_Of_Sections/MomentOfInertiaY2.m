function I = MomentOfInertiaY2(IVars, IMethod)
% function: MomentOfInertiaY2
%
% A function for calculating the second moment of inertia of a beam about 
% the y2 axis.
%
% Input:
%
%   vars - (n x p) column vector of values necessary for calculating the
%          second moment of inertia according to the method described in 
%          'method'. See the appropriate function for details of the
%          required variables, and exact format of vars.
%
%   method - string describing the method by which the second moment of 
%            inertia is to be calculated. These should correspond to the 
%            appropriate table in Roark's Formulas for Stress & Strain. If
%            not one of the listed functions, MomentOfArea will interpret
%            this as a function name and attmpt to evaluate the named
%            function using Ivars as the input parameters.
% 
% Output:
%
%   I - (n x 1) column vector of values of I, the second moment of area of 
%       the cross-section about the axis y2, calulated according to 'method'
%
    switch IMethod
        
        case '1.2'
            % Rectangular
            I = Table1r2I_2(IVars);
            
        case '1.3'
            % Hollow Rectangular
            I = Table1r3I_2(IVars);

        case '1.6'
            % I-Beam
            I = Table1r6I_2(IVars);

        case '1.14'
            % Solid Circle: same as I_1
            I = Table1r14I_2(IVars);

        case '1.15'
            % Hollow Circle: same as I_1
            I = Table1r15I_2(IVars);

        case '1.16'
            % Hollow Circle with very thin annulus: same
            % as I-1
            I = Table1r16I_2(IVars);

        otherwise

            I = feval(IMethod, IVars);
            
    end

end