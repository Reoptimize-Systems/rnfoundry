function span2val = Table1r3span2(vars)
% Table1r3span2: Calculates the maximum span of a beam with a
% hollow rectangular cross-section, in the axis 2
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
%   span2val - (n x 1) column vector of values of the maximum span of the
%   section in axis 2
%

    d = vars(:,1);

    span2val = d;
    
end