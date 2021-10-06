function [thetaT, omegaT] = prescribedmotomegatheta(t, simoptions)
% evaluates piecewise polynomials fitted to a time series of angular
% position and velocity values
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
%       pp_thetaT - a piecewice polynomials fitted to the angular positions
%                   to be interpolated
%       pp_omegaT - a piecewice polynomials fitted to the angular
%                   velocities to be interpolated
%
% Output
%
%   thetaT - the angular positions at the times in 't'
%
%   omegaT - the angular velocities at the times in 't'
%
%
% See also: prescribedmotodetorquefcn_ROTARY
%

% Created by Richard Crozier 2012-2015

    thetaT = ppuval(t, simoptions.pp_thetaT);
    
    omegaT = ppuval(t, simoptions.pp_omegaT);

end