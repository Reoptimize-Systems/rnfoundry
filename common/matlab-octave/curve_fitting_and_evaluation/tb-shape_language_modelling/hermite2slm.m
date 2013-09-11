function slm = hermite2slm(harray)
% hermite2slm: converts an array of coefficients in Hermite form to a slm
% usage: slm = hermite2slm(harray)
%
% arguments: (input)
%  harray - nx2 or nx3 array of coefficients
%           harray(:,1) == knots
%           harray(:,2) == function values
%
%           and if they are present:
%           harray(:,3) == first derivatives
%
%           If additional columns, they define the
%           higher order derivatives of a (2*k-3)
%           degree hermite interpolant.
%
% arguments: (output)
%  slm - slm struct, usable by slmeval and plotslm
%
%
% See also: slmset, slmengine, slmeval, ppval, slmfit, plotslm
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 2/6/07

% ------------------           licence           --------------------
% Copyright (c) 2007, John D'Errico
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
% ------------------           licence           --------------------

[n,k] = size(harray);

% create the basic struct
slm.form = 'slm';
slm.knots = harray(:,1);
% we don't know how it was created
slm.prescription = [];
% nor what it was derived from
slm.x = [];
slm.y = [];

if k == 2
  % linear Hermite
  slm.degree = 1;
  slm.coef = harray(:,2);
  
elseif k == 3
  % cubic Hermite
  slm.degree = 3;
  slm.coef = harray(:,2:3);
  
elseif k > 3
  % quintic Hermite or higher
  slm.degree = 2*k-3;
  slm.coef = harray(:,2:end);

else
  error 'harray is not in my Hermite standard form'
end

