function alphabeta = abc2alphabeta (abc)
% performs the clarke transformation on a set of three-phase quantities
%
% Syntax
%
% dqo = abc2alphabeta (abc)
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
%  abc - vector of three values, the quantities (e.g. currents or voltages)
%   in a balanced three phase system which are to be projected onto the
%   direct quadrature zero axes
%
% Output
%
% alphabeta - vector of three values containing the transformed quantities
%   on the rotating reference frame
%
% See also: alphabeta2abc.m, alphabeta2dq0.m, dq02alphabeta.m
%           abc2dq0.m, dq02abc.m
%


    if ~(isvector(abc) && numel(abc) == 3 && all(isreal(abc)))
        error('abc2alphabeta:badabcinput', 'abc must be a vector of three real values.')
    end
    
    rt3over2 = sqrt(3) / 2;
    
    alphabeta = (2/3) * [ 1,    -0.5,       -0.5;
                          0,    rt3over2,   -rt3over2;
                          0.5,  0.5,        0.5 ] ...
                      * abc(:);

end