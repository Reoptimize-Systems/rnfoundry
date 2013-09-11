function s = polyn2sym(polyn)
% polyn2sympoly: convert a regression polynomial from polyfitn into its symbolic toolbox form
% usage: sp = polyn2sym(polyn)
%
% arguments: (input)
%  polyn - a structure as returned from polyfitn
%
% arguments: (output)
%  s - A symbolic toolbox object
%
% After conversion into a symbolic toolbox form, any symbolic operations are
% now possible on the polynomial.
% 
% Requirement: The symbolic toolbox, as supplied by the MathWorks.
% http://www.mathworks.com/products/symbolic/functionlist.html
%
% See also: polyvaln, polyfit, polyval, sym
%
% Author: John D'Errico
% Release: 3.0
% Release date: 8/23/06

% ----------------          licence         -----------------
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
% ----------------          licence         -----------------

if exist('sym','file')~=2
  error 'Please obtain the symbolic toolbox from the MathWorks'
end

% initialize the returned argument as symbolic
s = sym(0);

% Unpack the fields of polyn for use
Varlist = polyn.VarNames;
Expon = polyn.ModelTerms;
Coef = polyn.Coefficients;
% how many terms?
nterms = length(Coef);

% Was the list of variable names empty?
% If so, then generate a list of names of my own.
nvars = size(polyn.ModelTerms,2);
if isempty(Varlist)
  Varlist={};
  for i = 1:nvars
    Varlist{i} = ['X',num2str(i)];
  end
end

% make the vars symbolic
for i = 1:nvars
  Varlist{i} = sym(Varlist{i});
end

% build the polynomial
for i = 1:nterms
  term = sym(Coef(i));
  for j = 1:nvars
    term = term*Varlist{j}^Expon(i,j);
  end
  
  % accumulate into s
  s = s+term;
end


