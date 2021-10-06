function [str, pfs] = sipre(val,sgf,pfx,trz)
% Convert a scalar numeric into an SI prefixed string. (metric/engineering)
%
% Syntax:
%  str = sipre(val)             % Four significant figures and prefix symbol.
%  str = sipre(val,sgf)         % Select significant figures, prefix symbol.
%  str = sipre(val,sgf,pfx)     % Select sig-figs, choose prefix symbol or name.
%  str = sipre(val,sgf,pfx,trz) % Select if decimal trailing zeros are required.
%
% Description
%
% Convert a scalar numeric value into a string. The value is shown in the
% string as a coefficient and an SI unit prefix, optimally chosen for
% readability. If the rounded |val|<10^-24 or |val|>=10^27 then E-notation
% is used, without a prefix.
%
% The following SI prefix string and symbols are used:
%
% Order  |1000^1 |1000^2 |1000^3 |1000^4 |1000^5 |1000^6 |1000^7 |1000^8 |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Name   | kilo  | mega  | giga  | tera  | peta  | exa   | zetta | yotta |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Symbol*|   k   |   M   |   G   |   T   |   P   |   E   |   Z   |   Y   |
%
% Order  |1000^-1|1000^-2|1000^-3|1000^-4|1000^-5|1000^-6|1000^-7|1000^-8|
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Name   | milli | micro | nano  | pico  | femto | atto  | zepto | yocto |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Symbol*|   m   |   u   |   n   |   p   |   f   |   a   |   z   |   y   |
%
%
% Input
%
%  val - NumericScalar, the value to be converted to string <str>.
%
%  sgf - (optional) NumericScalar, the significant figures in the
%    coefficient, default is 4.
%
%  pfx - (optional) LogicalScalar, true/false flag selecting if SI prefix
%    is a full name or symbol (e.g. 'mega' or 'M'). If true, the full name
%    is used, otherwise the symbol. Default is false.
%
%  trz -(optional)  LogicalScalar, true/false selecting if decimal trailing
%    zeros are required.
%
% Output
%
%  str - Input <val> as a string: coefficient + space character + SI prefix.
%
% str = sipre(val,*sgf,*pfx,*trz)
%
% See also SINUM BIPRE BINUM NUM2STR STR2NUM MAT2STR SSCANF SPRINTF 
%          ROUND60063 ROUND2SF ROUND2DP NUM2WORDS
%
% Examples
%
% sipre(10000)  OR  sipre(1e4)
%   ans = '10 k'
% sipre(10000,4,true)
%   ans = '10 kilo'
% sipre(10000,4,false,true)
%   ans = '10.00 k'
%
% ['Power: ',sipre(200*1000^2,2,true),'watt']
%   ans = 'Power: 200 megawatt'
%
% sipre(-5.555e9,2) % Rounds significant figures correctly.
%   ans = '-5.6 G'
%
% sprintf('Clock frequency is %shertz.',sipre(1234567890,5,true))
%   ans = 'Clock frequency is 1.2346 gigahertz.'
%
% sipre(sinum('9 T'))
%   ans = '9 T'
%



% Copyright (c) 2014, Stephen Cobeldick
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
%     * Neither the name of the  nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
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


% ### Input Wrangling ###
%
if nargin<4
    trz = false;
else
    assert(islogical(trz)&&isscalar(trz),'Fourth input <trz> must be a logical scalar.')
end
if nargin<3
    pfx = false;
else
    assert(islogical(pfx)&&isscalar(pfx),'Third input <pfx> must be a logical scalar.')
end
if nargin<2
    sgf = 4;
else
    assert(isnumeric(sgf)&&isscalar(sgf),'Second input <sgf> must be a numeric scalar.')
    sgf = double(uint8(sgf));
end
assert(isnumeric(val)&&isscalar(val)&&isreal(val),'First input <val> must be a real numeric scalar.')
val = double(val);
%
if trz && sgf>1
    fmt = '%#.*g %s';
else
    fmt = '%.*g %s';
end
%
% ### Generate String ###
%
if isfinite(val)
    % Calculate coefficient value:
    xpt = rem(min(9,max(-9,[0;1]+floor(log10(abs(val))/3))),9);
    cof = val.*1000.^-xpt;
    % Round coefficient value:
    ord = 1+floor(log10(abs(cof)));
    if val~=0
        cof = 10.^(ord-sgf).*round(cof.*10.^(sgf-ord));
    end
    % Select prefix symbol/name:
    pfc = {'yocto','zepto','atto','femto','pico','nano','micro','milli',...
        '','kilo', 'mega', 'giga','tera', 'peta','exa', 'zetta','yotta';...
           'y',    'z',    'a',   'f',    'p',   'n',   'u',    'm',...
        '','k',    'M',    'G',   'T',    'P',   'E',   'Z',    'Y'};
    idx = 1+any(abs(cof)==[1000;1]);
    pfs = pfc{2-pfx,9+xpt(idx)};
    % Convert to string (without prefix || digits>whole part):
    if abs(ord(idx))>4 || floor(log10(abs(cof(idx)))-sgf)<-1
        str = sprintf(fmt,sgf,cof(idx),pfs);
    else % (digits<=whole part)
        str = sprintf('%.0f %s',cof(idx),pfs);
    end
else
    str = sprintf('%f ',val);
end
%
end
%----------------------------------------------------------------------END:sipre