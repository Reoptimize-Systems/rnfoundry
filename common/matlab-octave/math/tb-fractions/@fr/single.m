function S = single(FRAC)
% fr/single: conversion from a fraction object to a single.
% usage: S = single(FRAC);
% 
% arguments: (input)
%  FRAC - a fraction object
%
% arguments: (output)
%  S  - the single precision representation of FRAC.
% 
% Example:
%  m = fr(1,3);
%  single(m)
%  ans = 
%     0.3333
%
%
%  See also: double, fr

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/single as a template)

S=zeros(size(FRAC));
S(:)=single([FRAC(:).whole])+single([FRAC(:).numer])./single([FRAC(:).denom]);

end
