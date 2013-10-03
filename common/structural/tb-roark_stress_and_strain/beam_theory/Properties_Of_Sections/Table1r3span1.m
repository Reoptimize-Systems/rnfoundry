function span1val = Table1r3span1(vars)
% Table1r3span1: Calculates the maximum span of a beam with a
% hollow rectangular cross-section, in the axis 1
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
%   span1val - (n x 1) column vector of values of the maximum span of the
%   section in axis 1
%

    b = vars(:,2);

    span1val = b;
    
end