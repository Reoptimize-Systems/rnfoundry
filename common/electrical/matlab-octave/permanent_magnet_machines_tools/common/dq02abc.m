function iabc = dq02abc (idqo, theta)
% performs the inverse direct quadrature zero transformation on a set of
% d-q axis currents
%
% Syntax
%
% iabc = invdqo (idqo, theta)
%
% Description
%
% The direct quadrature zero (or dq0 or abc2dq0) transformation or
% zero direct quadrature (or 0dq or odq) transformation is a mathematical
% transformation used to simplify the analysis of three-phase circuits.
% This function performs the inverse abc2dq0 transform to obtain the equivalent
% three-phase currents. See 'abc2dq0' for more information.
%
% Input
%
% idqo - vector of two or three values, the d-q currents, or d-q-o
%   currents. If only two values are supplied, these are assumed to be the
%   d-q currents, and the third value is set to zero.
%
% theta - angle between the d axis and first three phase current
%
% Output
%
% iabc - vector of three values containing the three phase currents
%   obtained from the inverse transform of the d-q axis values
%

    if ~(isvector(idqo) && all(isreal(idqo)))
        error('DQO:badcurrentinput', 'idqo must be a real valued vector.')
    end
    
    if numel(idq) == 2
        idqo(3) = 0;
    elseif numel(idq) ~= 3
        error('DQO:badcurrentinput', 'idqo must be a vector of 2 or three values')
    end
    
    rt2over2 = sqrt(2) / 2;
    
    iabc = sqrt(2/3) * [ cos(theta),          sin(theta),          rt2over2;
                         cos(theta - 2*pi/3), sin(theta - 2*pi/3), rt2over2;
                         cos(theta + 2*pi/3), sin(theta + 2*pi/3), rt2over2] * idqo(:);
                     
end