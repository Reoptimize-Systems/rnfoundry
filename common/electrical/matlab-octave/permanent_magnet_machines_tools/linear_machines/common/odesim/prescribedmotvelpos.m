function [xT, vT] = prescribedmotvelpos(t, simoptions)
% evaluates piecewise polynomials fitted to a time series of position and
% velocity values
%
% Syntax
%
% [thetaT, omegaT] = prescribedmotomegatheta(t, simoptions)
%
% Input
%
%   t - the time points at which to evaluate the polynomials
%
%   simoptions - a structure containng at least the fields
%
%       pp_xT - a piecewice polynomials fitted to the positions to be
%               interpolated
%       pp_vT - a piecewice polynomials fitted to the velocities to be
%               interpolated
%
% Output
%
%   xT - the positions at the times in 't'
%
%   vT - the velocities at the times in 't'
%
%
% See also: prescribedmotode_linear, prescribedmotodeforcefcn_linear
%

% Created by Richard Crozier 2013


    xT = ppval(simoptions.pp_xT, t);
    
    vT = ppval(simoptions.pp_vT, t);

end