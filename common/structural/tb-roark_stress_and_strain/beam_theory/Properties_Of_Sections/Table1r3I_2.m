function I = Table1r3I_2(IVars)
% Table1r3I_1: Calculates the second moment of inertia of a beam 
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
%   I - (n x 1) column vector of values of I, the moment of inertia for a
%       beam with a hollow rectangular cross-section
%

    if size(IVars,2) == 4
        
        d = IVars(:,1);
        b = IVars(:,2);
        di = IVars(:,3);
        bi = IVars(:,4);

        I = ((b.^3) .* d - (bi.^3).*di) ./ 12;
    else
        error('Dimensions of Input variables incorrect. Table1r2I_2 accepts a (n x 2) column vector of values')
    end
     
end