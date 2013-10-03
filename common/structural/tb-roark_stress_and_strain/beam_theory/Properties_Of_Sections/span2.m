function span2val = span2(vars, Method)
% function: span2: maximum span of a section in the 2nd axis (axis 2 as
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
%   span2val - (n x 1) column vector of values of span2val, the maximum
%   span of the section in axis 2, calulated according to 'method'
%

    switch Method
        
        case '1.2'
            % Rectangular
            span2val = Table1r2span2(vars);
            
        case '1.3'
            % Hollow Rectangle
            span2val = Table1r3span2(vars);
        
        case '1.6'
            % I-Beam
            span2val = Table1r6span2(vars);

        case '1.14'
            % Solid Circle
            span2val = Table1r14span2(vars);

        case '1.15'
            % Hollow Circle
            span2val = Table1r15span2(vars);

        case '1.16'
            % Hollow Circle with very thin annulus
            span2val = Table1r16span2(vars);

        otherwise

            span2val = feval(Method, vars);
            
    end
end