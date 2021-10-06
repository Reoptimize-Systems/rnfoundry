function span1val = span1(vars, Method)
% function: span1: maximum span of a section in the 1st axis (axis 1 as
% defined in Roark's Formulas for Stress and Strain 6th Edition)
%
% Input:
%
%   vars - (n x p) column vector of values necessary for calculating the
%          maximum span according to the method described in 
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
%   span1val - (n x 1) column vector of values of span1val, the maximum
%   span of the section in axis 1, calulated according to 'method'
%

    switch Method
        
        case '1.2'
            % Rectangular
            span1val = Table1r2span1(vars);
            
        case '1.3'
            % Hollow Rectangle
            span1val = Table1r3span1(vars);
        
        case '1.6'
            % I-Beam
            span1val = Table1r6span1(vars);

        case '1.14'
            % Solid Circle
            span1val = Table1r14span1(vars);

        case '1.15'
            % Hollow Circle
            span1val = Table1r15span1(vars);

        case '1.16'
            % Hollow Circle with very thin annulus
            span1val = Table1r16span1(vars);

        otherwise

            span1val = feval(Method, vars);
            
    end
end