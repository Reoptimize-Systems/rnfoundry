function dqo = abc2dq0 (abc, theta)
% performs the clarke transformation on a set of three-phase quantities
%
% Syntax
%
% dqo = abc2dq0 (abc, theta)
%
% Description
%
% The direct quadrature zero (or dq0 or dqo) transformation or zero direct
% quadrature (or 0dq or odq) transformation is a mathematical
% transformation used to simplify the analysis of three-phase circuits. In
% the case of balanced three-phase circuits, application of the dq0
% transform reduces the three AC quantities to two DC quantities.
% Simplified calculations can then be carried out on these imaginary DC
% quantities before performing the inverse transform to recover the actual
% three-phase AC results. It is often used in order to simplify the
% analysis of three-phase synchronous machines or to simplify calculations
% for the control of three-phase inverters. The dq0 transform is similar to
% the transform first proposed in 1929 by R.H. Park. In fact, the dq0
% transform is often referred to as Park s transformation, although there
% are in fact important differences.
%
% Input
%
% abc - vector of three values, the quantities (e.g. currents or voltages)
%   in a balanced three phase system which are to be projected onto the
%   direct quadrature zero axes
%
% theta - angle between the d axis and first current in abc, this is given
%   by omega*t in a time varying system 
%
% Output
%
% dqo - vector of three values containing the transformed quantities on the
%   direct quadrature zero axes.
%

    if ~(isvector(abc) && numel(abc) == 3 && all(isreal(abc)))
        error('DQO:badabcinput', 'abc must be a vector of three real values.')
    end
    
    rt2over2 = sqrt(2) / 2;
    
    dqo = sqrt(2/3) * [ cos(theta), cos(theta - 2*pi/3), cos(theta + 2*pi/3);
                        sin(theta), sin(theta - 2*pi/3), sin(theta + 2*pi/3);
                        rt2over2,   rt2over2,            rt2over2] ...
                    * abc(:);

end