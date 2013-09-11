%{

Author: John D'Errico
Release: 1.2
Release date: 2/20/06\cf0 \

----------------          licence         -----------------
Copyright (c) 2007, John D'Errico
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright 
      notice, this list of conditions and the following disclaimer in 
      the documentation and/or other materials provided with the distribution
      
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.
----------------          licence         -----------------

What follows are example usages of polyfitn, polyvaln, and poly2sympoly.

%}

% Fit a 1-d model to cos(x). We only need the even order terms.
x = -2:.1:2;
y = cos(x);
p = polyfitn(x,y,'constant x^2 x^4 x^6')

% Conversion to a sympoly. If nothing else, its a nice way to display the model.
polyn2sympoly(p)

% Evaluate the regression model at some set of points
polyvaln(p,[0 .5 1])

%%

% A surface model in 2-d, with all terms up to third order.
% Use lots of data.
n = 1000;
x = rand(n,2);
y = exp(sum(x,2)) + randn(n,1)/100;
p = polyfitn(x,y,3)
polyn2sympoly(p)

% Evaluate on a grid and plot:
[xg,yg]=meshgrid(0:.05:1);
zg = polyvaln(p,[xg(:),yg(:)]);
surf(xg,yg,reshape(zg,size(xg)))
hold on
plot3(x(:,1),x(:,2),y,'o')
hold off

%%

% A linear model, but with no constant term, in 2-d
uv = rand(100,2);
w = sin(sum(uv,2));
p = polyfitn(uv,w,'u, v');
polyn2sympoly(p)

%%

% A model with various exponents.
% Also with only 1 variable, x & y may be row or column vectors.
x = 1:10;
y = 3 + 2./x + sqrt(x) + randn(size(x))/100;
p = polyfitn(x,y,'constant x^-1 x^0.5');
polyn2sympoly(p)

xi = 1:.1:10;
yi = polyvaln(p,xi);
plot(x,y,'ro',xi,yi,'b-')


