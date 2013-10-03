function A = Table1r3Area(vars)
% Table1r3Area: Calculates the cross-sectional area of a beam 
% with a hollow rectangular cross-section, as calculated in 'Roark's 
% Formulas Stress & Strain 6th edition' in table 1, page 62 row 3.
%
%
%                  2
%                  .
%            ______.______
%        ^  |  _________  |
%    ^   ¦  | |    .    | |
%    ¦   ¦  | |    .    | |
% di ¦d  ¦  |.|.........|.|......... 1
%    ¦   ¦  | |    .    | |
%    ¦   ¦  | |    .    | |
%    v   ¦  | |_________| |
%        v  |______.______|
%           <------.------>
%               b  .
%                  .
%              <------>
%                 bi
%
% Input: 
%
%   IVars - (n x 4) matrix of values as below
%           IVars(:,1) - d,rectangle height
%           IVars(:,2) - b, width
%           IVars(:,3) - di, inner rectangle void height
%           IVars(:,4) - bi, inner rectangular void width
%
% Output:
%
%   A - (n x 1) column vector of values of A, the cross-sectional area of a
%       beam with a hollow rectangular cross-section
%

    d = vars(:,1);
    b = vars(:,2);
    di = vars(:,3);
    bi = vars(:,4);
    
    A = (b .* d) - (bi .* di);

end