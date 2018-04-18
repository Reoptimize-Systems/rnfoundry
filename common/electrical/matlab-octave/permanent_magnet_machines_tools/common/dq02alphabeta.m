function alphabeta = dq02alphabeta (dq0, theta)
% performs the inverse dq0 (Park) transformation on a set of values
%
% Syntax
%
% alphabeta = dq02alphabeta (dqo, theta)
%
% Description
%
% modifies the voltages in d,q rotating reference frame in a two phase
% orthogonal system to get the alpha-beta values.
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
%  alphabeta - vector of three values containing the three phase quantities
%   obtained from the inverse transform of the d-q axis values
%
% See also: alphabeta2dq0.m, abc2alphabeta.m, alphabeta2abc.m
%           abc2dq0.m, dq02abc.m
%

    if ~(isvector(dq0) && all(isreal(dq0)))
        error('dq02alphabeta:baddqoinput', 'dqo must be a real valued vector of length 2 or 3.')
    end
    
    alphabeta = [ dq0(1) .* cos(theta) - dq0(2) .* sin(theta);
                  dq0(1) .* sin(theta) + dq0(2) .* cos(theta); ];
                     
end