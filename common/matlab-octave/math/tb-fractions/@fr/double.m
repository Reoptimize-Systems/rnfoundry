function D = double(FRAC)
% fr/double: conversion from a fraction object to a double.
% usage: D = double(FRAC);
% 
% arguments: (input)
%  FRAC - a fraction object
%
% arguments: (output)
%  D  - the double precision approximation of FRAC.
% 
% Example:
%  m = fr(1,3);
%  double(m)
%  ans = 
%     0.3333
%
%
%  See also: single, fr

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/double as a template)

D=zeros(size(FRAC));
D(:)=double([FRAC(:).whole])+double([FRAC(:).numer])./double([FRAC(:).denom]);

end
