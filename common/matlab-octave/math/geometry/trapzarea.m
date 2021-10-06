function area = trapzarea(a, b, h)
% calculates the area of a trapesium
%
% Syntax
%
% area = trapzarea(a, b, h)
%
% Description 
%
%        <-- a -->
%        _________
%       /          \       ^
%      /            \      : h
%     /              \     :
%    /________________\    v
%    <------ b ------>
%
%
    area = h .* ((a + b) ./ 2);

end