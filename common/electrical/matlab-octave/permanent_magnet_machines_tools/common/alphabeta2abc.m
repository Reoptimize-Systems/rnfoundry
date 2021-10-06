function abc = alphabeta2abc (alphabetagamma)
% performs the inverse clarke transformation on a set of three-phase quantities
%
% Syntax
%
% abc = alphabeta2abc (alphabetagamma)
%
% Description
%
% The alpha-beta transformation (also known as the Clarke transformation)
% is a mathematical transformation employed to simplify the analysis of
% three-phase circuits. Conceptually it is similar to the dq0
% transformation.
%
% Input
%
%  dqo - vector of two or three values, the alpha-beta values, or
%   alpha-beta-gamma values. If only two values are supplied, these are
%   assumed to be the alpha-beta values, and the third value is set to
%   zero.
%
% Output
%
% alphabetagamma - vector of three values containing the three phase
%  quantities obtained from the inverse transform of the alpha-beta axis
%  values
%
% See also: abc2alphabeta.m, alphabeta2dq0.m, dq02alphabeta.m
%           abc2dq0.m, dq02abc.m
%


    if ~(isvector(alphabetagamma) && all(isreal(alphabetagamma)))
        error('alphabetagamma:baddqoinput', 'alphabetagamma must be a real valued vector.')
    end
    
    if numel(alphabetagamma) == 2
        alphabetagamma(3) = 0;
    elseif numel(alphabetagamma) ~= 3
        error('alphabetagamma:baddqoinput', 'alphabetagamma must be a vector of 2 or three values')
    end
    
    rt3over2 = sqrt(3) / 2;
    
    abc = [ 1,     0,          1;
            -0.5,  rt3over2,   1;
            -0.5,  -rt3over2,  1 ] * alphabetagamma(:);

end