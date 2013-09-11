function sp = polyn2sympoly(polyn)
% polyn2sympoly: convert a regression polynomial from polyfitn into a sympoly
% usage: sp = polyn2sympoly(polyn)
%
% arguments: (input)
%  polyn - a structure as returned from polyfitn
%
% arguments: (output)
%  sp - A sympoly object
%
% After conversion into a sympoly, any symbolic operations are
% now possible on this form.
% 
% Requirement: The sympoly toolbox, as found on Matlab Central.
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=9577&objectType=FILE
%
% See also: polyvaln, polyfit, polyval, sympoly
%
% Author: John D'Errico
% Release: 1.0
% Release date: 2/19/06

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

if exist('sympoly','file')~=2
  error 'Please download the sympoly toolbox from Matlab Central'
end

% Copy over the fields of polyn into sp
sp.Var = polyn.VarNames;
sp.Exponent = polyn.ModelTerms;
sp.Coefficient = polyn.Coefficients(:);

% Was the list of variable names empty?
% If so, then generate a list of names of my own.
if isempty(sp.Var)
  p = size(polyn.ModelTerms,2);
  varlist={};
  for i = 1:p
    varlist{i} = ['X',num2str(i)];
  end
  sp.Var = varlist;
end

% turn the struct into a sympoly
sp = sympoly(sp);


