function span2val = Table1r6span2(vars)
% Table1r6span2: Calculates the maximum span in axis 2 of a beam with a I-beam
% cross-section, (wide-flange beam with equal flanges) as calculated in
% 'Roark's Formulas Stress & Strain 6th edition' in table 1, page 64 row 6. 
%
% Input: 
%   
%   vars - (n x 4) matrix of values as below
%           vars(:,1) - b, width of the I-beam flanges
%           vars(:,2) - t, the thickness of the flanges
%           vars(:,3) - tw, the thickness of the I-Beam vertical part
%           vars(:,4) - d, height of the I-Beam vertical part (not
%           including flange thickness)
%
%   outer radius, the second contains values of t, the annulus thickness.
%
% Output:
%
%   span2val - (n x 1) column vector of values of the maximum span of the
%   section in axis 1
%

    t = vars(:,2);
    d = vars(:,1);
    
    span2val = d + t;

end