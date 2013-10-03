function I = Table1r6I_1(IVars)
% function: Table1r6I_1
% 
% Calculates the moment of inertia of a beam with a I-beam
% cross-section, (wide-flange beam with equal flanges) as calculated in
% 'Roark's Formulas Stress & Strain 6th edition' in table 1, page 64 row 6. 
%
% Input: 
%   
%   IVars - (n x 4) matrix of values as below
%           IVars(:,1) - b, width of the I-beam flanges
%           IVars(:,2) - t, the thickness of the flanges
%           IVars(:,3) - tw, the thickness of the I-Beam vertical part
%           IVars(:,4) - d, height of the I-Beam vertical part (not
%           including flange thickness)
%
%   outer radius, the second contains values of t, the annulus thickness.
%
% Output:
%
%   I - (n x 1) column vector of values of I, the moment of inertia for a
%   beam with an I-beam cross-section
%
    if size(IVars,2) == 4
        
        b = IVars(:,1);
        t = IVars(:,2);
        tw = IVars(:,3);
        d = IVars(:,4);

        I = (b.*(d + 2.*t).^3) ./ 12 - (b - tw) .* d.^3 ./ 12;

    else
        error('Dimensions of Input variables incorrect. Table1r6I_2 accepts a (n x 4) column vector of values')
    end

end