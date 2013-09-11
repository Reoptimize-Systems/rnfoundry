function idqo = dqo(iabc, theta)
% performs the direct–quadrature–zero transformation on a set of
% three-phase currents
%
% Syntax
%
% idqo = dqo(iabc, theta)
%
% Description
%
% The direct–quadrature–zero (or dq0 or dqo) transformation or
% zero–direct–quadrature (or 0dq or odq) transformation is a mathematical
% transformation used to simplify the analysis of three-phase circuits. In
% the case of balanced three-phase circuits, application of the dqo
% transform reduces the three AC quantities to two DC quantities.
% Simplified calculations can then be carried out on these imaginary DC
% quantities before performing the inverse transform to recover the actual
% three-phase AC results. It is often used in order to simplify the
% analysis of three-phase synchronous machines or to simplify calculations
% for the control of three-phase inverters. The dqo transform is similar to
% the transform first proposed in 1929 by R.H. Park. In fact, the dqo
% transform is often referred to as Park’s transformation, although there
% are in fact important differences.
%
% Input
%
% iabc - vector of three values, the currents in a balanced three phase
%   system which are to be projected onto the direct–quadrature–zero axes
%
% theta - angle between the d axis and first current in iabc
%
% Output
%
% idqo - vector of three values containing the transformed currents on the
%   direct–quadrature–zero axes. 

    if ~(isvector(iabc) && numel(iabc) == 3 && all(isreal(iabc)))
        error('DQO:badcurrentinput', 'i must be a vector of three real values.')
    end
    
    rt2over2 = sqrt(2) / 2;
    
    idqo = sqrt(2/3) * [cos(theta), cos(theta - 2*pi/3), cos(theta + 2*pi/3);
                        sin(theta), sin(theta - 2*pi/3), sin(theta + 2*pi/3);
                        rt2over2,   rt2over2,            rt2over2] * iabc(:);

end