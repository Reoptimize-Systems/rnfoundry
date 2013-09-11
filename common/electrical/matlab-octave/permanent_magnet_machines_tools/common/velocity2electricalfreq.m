function fe = velocity2electricalfreq(v, polewidth)
% converts a linea machine velocity to its electrical frequency at this
% speed
%
% Syntax
%
% fe = velocity2electricalfreq(v, polewidth)
%
% Input
%
%   v - the relative velocity of the two machine parts
%
%   polewidth - the width of one electrical pole of the machine
%
% Output
%
%   fe - electrical frequency at the specified velocity
%

    % s = d/t
    % st = d
    % t = d / s
    % f = s / d
    fe = abs(v) ./ (2 * polewidth);

end
