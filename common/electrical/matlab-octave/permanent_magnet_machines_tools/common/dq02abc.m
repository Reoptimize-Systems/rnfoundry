function abc = dq02abc (dqo, theta, startaxisd)
% performs the inverse dq0 (Park) transformation on a set of values
%
% Syntax
%
% abc = invdqo (dqo, theta)
%
% Description
%
% The direct quadrature zero (or dq0 or abc2dq0) transformation or zero
% direct quadrature (or 0dq or odq) transformation is a mathematical
% transformation used to simplify the analysis of three-phase circuits.
% This function performs the inverse abc2dq0 transform to obtain the
% equivalent three-phase currents. See 'abc2dq0' for more information.
%
% Input
%
%  dqo - vector of two or three values, the d-q values, or d-q-o
%   values. If only two values are supplied, these are assumed to be the
%   d-q values, and the third value is set to zero.
%
%  theta - angle between the d axis and first three phase quantity
%
% Output
%
%  abc - vector of three values containing the three phase quantities
%   obtained from the inverse transform of the d-q axis values
%

    if ~(isvector(dqo) && all(isreal(dqo)))
        error('dq02abc:baddqoinput', 'dqo must be a real valued vector.')
    end
    
    if numel(dqo) == 2
        dqo(3) = 0;
    elseif numel(dqo) ~= 3
        error('dq02abc:baddqoinput', 'dqo must be a vector of two or three values')
    end
    
    pi2v3 = 2*pi/3;
    
    if startaxisd
        
        abc = [ cos(theta),           -sin(theta),          1;
                cos(theta - pi2v3),   -sin(theta - pi2v3),  1;
                cos(theta + pi2v3),   -sin(theta + pi2v3),  1 ] ...
              * dqo(:);
          
    else
        
        abc = sqrt(2/3) * [ sin(theta),            cos(theta),          1;
                            sin(theta - pi2v3),    cos(theta - pi2v3),  1;
                            sin(theta + pi2v3),    cos(theta + pi2v3),  1 ] ...
                        * dqo(:);
                    
    end
    
end