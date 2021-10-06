function P = corelossfcn_AM(v, x, hxslm, cxslm, exslm)
% evaluates the soft ferrite losses in the core from partially
% precalculated loss functions
%
% Syntax
%
% P = corelossfcn(v, x, hxslm, cxslm, exslm)
%
% Input
%
%   v - relative velocity of the machine parts
%
%   x - relative position of the machine parts
%
%   hxslm - slm object fitted to the hysteresis loss function over two
%     Poles
%
%   cxslm - slm object fitted to the eddy current loss function over two
%     Poles
%
%   exslm - slm object fitted to the excess loss function over two Poles
%
%
% Output
%
%   P - the instantaneous power loss in the iron core
%
    
% Copyright Richard Crozier 2012

    P = abs(v .* periodicslmeval(x, hxslm, 0, false) ...
        + realpow(v,2) .* periodicslmeval(x, cxslm, 0, false) ...
        + realpow(sqrt(v),3) .* periodicslmeval(x, exslm, 0, false));
    
    % smooth the function around zero velocity for the ode solvers. Without
    % introducing a gradient below a threshold the forces calculated from
    % the losses flip from positive to negative about zero velocity. This
    % causes problems for the ode solvers which attempt to take smaller and
    % smaller time steps to maintain the acuracy of the solution.
    vthresh = 0.001;
    % the function v^2 / vthresh^2 is a quadratic which rises from zero to
    % 1 between v = 0 and v = vthresh. 
    P(v < vthresh) = P(v < vthresh) .* v(v < vthresh).^2 ./ (vthresh^2); 
    
end