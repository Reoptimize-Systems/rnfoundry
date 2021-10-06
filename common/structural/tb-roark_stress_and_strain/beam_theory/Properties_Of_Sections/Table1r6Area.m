function A = Table1r6Area(vars)
% Table1r6Area: Calculates the cross-sectional area of a beam with a I-beam
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
%   A - (n x 1) column vector of values of A, the cross-sectional area of a
%   beam with an I-beam cross-section
%

    b = vars(:,1);
    t = vars(:,2);
    tw = vars(:,1);
    d = vars(:,2);
    
    A = (2 .* b .* t) + (tw .* d);

end