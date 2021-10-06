function K = Table20r16K(vars)
% Table20r16K: Calculates the cross-sectional area of a beam with a hollow,
% thin walled rectangular cross-section, as calculated in 'Roark's Formulas
% Stress & Strain 6th edition' in table 20, page 352 row 16.
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
%   K - (n x 1) column vector of values of K, the torsion constant of a
%   beam with a hollow, thin walled, rectangular cross-section
%

    d = vars(:,1);
    b = vars(:,2);
    di = vars(:,3);
    bi = vars(:,4);
    
    t = (b - bi) ./ 2;
    
    t1 = (d - di) ./ 2;
    
    K = 2 .* t .* t1 .* (b - t).^2 .* (d - t1).^2 ./ (b .* t + d .* t1 - t.^2 - t1.^2);
    
end