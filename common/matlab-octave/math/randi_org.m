function R = randi_org(imax, n)
%RANDI_ORG Random (www.random.org) integers from a uniform discrete distribution.
%
%   R = RANDI_ORG returns the quota of available random numbers.
%   R = RANDI_ORG(IMAX,N) returns an N-by-N matrix containing random
%   integer values drawn from the discrete uniform distribution on 1:IMAX.
%   RANDI_ORG(IMAX,[M,N]) returns an M-by-N matrix.
%   RANDI_ORG(IMAX,[M,N,P,...]) returns an
%   M-by-N-by-P-by-... array.
%   RANDI_ORG(IMAX,SIZE(A)) returns an array the same size as A.
%
%   R = RANDI_ORG([IMIN,IMAX],...) returns an array containing integer
%   values drawn from the discrete uniform distribution on IMIN:IMAX.
%
%   Note: The size inputs M, N, P, ... should be positive integers.
%
%   Note: the function connects via http to www.random.org. It depends on
%         a working internet connection.
%
%   Examples:
%
%      Generate integer values from the uniform distribution on the set 1:10.
%         r = randi_org(10,[100,1]);
%
%   See also RAND, RANDN, RANDSTREAM, RANDSTREAM/RANDI, RANDSTREAM.GETDEFAULTSTREAM.

%   Copyright (c) 2010, Giampiero Salvi
%   All rights reserved.
%
%   Redistribution and use in source and binary forms, with or without 
%   modification, are permitted provided that the following conditions are 
%   met:
%
%       * Redistributions of source code must retain the above copyright 
%         notice, this list of conditions and the following disclaimer.
%       * Redistributions in binary form must reproduce the above copyright 
%         notice, this list of conditions and the following disclaimer in 
%         the documentation and/or other materials provided with the distribution
%       
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%   POSSIBILITY OF SUCH DAMAGE.

if nargin == 0
    [result, status] = urlread('http://www.random.org/quota/?format=plain');
    R = str2num(result);
    return
end

if isscalar(imax)
    imin = 0;
else
    if length(imax(:))>2
        error('imax can at most hold two values');
    else
        imin = imax(1);
        imax = imax(2);
    end
end

if any(n<=0)
    error('randi_org:dimOutOfRange', 'n should contain positive integers');
end

if isscalar(n)
    sizes = [n n];
else
    sizes = n;
end

tot = prod(sizes(sizes>0));

% www.random.org restricts queries to 10000 numbers
maxlen = 10000;
nqueries = ceil(tot/maxlen);
R = zeros(tot,1);
for qu = 1:nqueries
    len = maxlen;
    if qu==nqueries
        len = rem(tot, maxlen);
    end
    places = (qu-1)*maxlen + (1:len);
    [result, status] = urlread(['http://www.random.org/integers/?num=' num2str(len) '&min=' num2str(imin) '&max=' num2str(imax) '&col=1&base=10&format=plain&rnd=new']);
    if status~=1
        error('randi_org:queryFail', 'command failed with code %d', status);
    end
    R(places) = str2num(result);
end
R = reshape(R, sizes);
