function A = CSArea(vars, Method)
% function: CSArea: calculates the cross-sectional area of a section as
% tabulated in Roark's Formulas for Stress and Strain 6th Edition
%
% Input:
%
%   vars - (n x p) column vector of values necessary for calculating the
%          section area according to the method described in 
%          'method'. See the appropriate function for details of the
%          required variables, and exact format of vars.
%
%   method - string describing the method by which the 
%            area is to be calculated. These should correspond to the 
%            appropriate table in Roark's Formulas for Stress & Strain. If
%            not one of the listed functions, CSArea will interpret
%            this as a function name and attmpt to evaluate the named
%            function using vars as the input parameters.
% 
% Output:
%
%   A - (n x 1) column vector of values of A, the cross-sectional area of
%       the section, calulated according to 'method'
%

    switch Method
        
        case '1.2'
            % Rectangular
            A = Table1r2Area(vars);
            
        case '1.3'
            % Hollow Rectangle
            A = Table1r3Area(vars);
        
        case '1.6'
            % I-Beam
            A = Table1r6Area(vars);

        case '1.14'
            % Solid Circle
            A = Table1r14Area(vars);

        case '1.15'
            % Hollow Circle
            A = Table1r15Area(vars);

        case '1.16'
            % Hollow Circle with very thin annulus
            A = Table1r16Area(vars);

        otherwise

            A = feval(Method, vars);
            
    end
end