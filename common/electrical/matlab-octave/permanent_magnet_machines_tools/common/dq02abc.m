function abc = dq02abc (dqo, theta)
% performs the inverse direct quadrature zero transformation on a set of quantities
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
% dqo - vector of two or three values, the d-q values, or d-q-o
%   values. If only two values are supplied, these are assumed to be the
%   d-q values, and the third value is set to zero.
%
% theta - angle between the d axis and first three phase quantity
%
% Output
%
% abc - vector of three values containing the three phase quantities
%   obtained from the inverse transform of the d-q axis values
%

    if ~(isvector(dqo) && all(isreal(dqo)))
        error('DQO:baddqoinput', 'dqo must be a real valued vector.')
    end
    
    if numel(idq) == 2
        dqo(3) = 0;
    elseif numel(idq) ~= 3
        error('DQO:baddqoinput', 'dqo must be a vector of 2 or three values')
    end
    
    rt2over2 = sqrt(2) / 2;
    
    abc = sqrt(2/3) * [ cos(theta),          sin(theta),          rt2over2;
                        cos(theta - 2*pi/3), sin(theta - 2*pi/3), rt2over2;
                        cos(theta + 2*pi/3), sin(theta + 2*pi/3), rt2over2 ] ...
                    * dqo(:);
                     
end