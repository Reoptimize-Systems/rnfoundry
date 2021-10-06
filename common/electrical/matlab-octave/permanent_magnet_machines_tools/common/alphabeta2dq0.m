function dq0 = alphabeta2dq0 (alphabeta, theta)
% performs the dq0 (Park) projection on a set of values
%
% Syntax
%
% alphabeta = dq02alphabeta (dqo, theta)
%
% Description
%
% Converts a set of alpha-beta pquntities representing a rotating
% three-phase system of values to the d-q representation.
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
% See also: dq02alphabeta.m, abc2alphabeta.m, alphabeta2abc.m,
%           abc2dq0.m, dq02abc.m
%

    if ~(isvector(alphabeta) && all(isreal(alphabeta)))
        error('alphabeta2dq0:baddqoinput', 'alphabeta must be a real valued vector of length 2.')
    end
    
    dq0 = [  alphabeta(1) .* cos(theta) + alphabeta(2) .* sin(theta);
            -alphabeta(1) .* sin(theta) + alphabeta(2) .* cos(theta);
            0 ];
                     
end