function length = rectcoilmtl(depth, height, thickness)
% calculates the mean turn length in a rectangular coil winding
%
% Syntax
% 
% length = rectcoilwirelen(depth, height, thickness)
% 
% Input
%                  depth
%         :<------------------->:
%        _:_____________________:_      :
%       /  _____________________  \ ....v..
%      |  |          :          |  |       height
%      |  |__________v__________|  |.......
%       \_________________________/     ^
%                    ^                  :
%                    : thickness
%
%   depth - the length of the internal rectangle
%
%   height - the height of the inner rectangle
%
%   thickness - the coil winding thickness
% 
% Output
%
%   length - the mean turn length

    % we will consider a rectangular coil a special case of an isosceles
    % trapezoidal coil winding, where the bottom and top parallel edges are
    % the same length
    length = isotrapzcoilmtl(depth, depth, height, thickness);
                                         
end