%% Fractions Toolbox
% This demo file shows how one might use the Fractions Toolbox
%
% The code was inspired by, and to a large extent based upon, John
% D'Errico's Variable Precision Integer Toolbox.
%
% Fractions are stored in the form K + N / D where N>=0 and D>=0.
%
% The toolbox was designed to be very flexible in the types of objects that
% can be used, for example doubles, integers and vpi numbers.  The fraction
% components can in principle be any objects provided certain minimal
% requirements are satisfied.
%
%
% Author: Ben Petschel  30/7/2009

%% The creator for fraction objects is fr
help fr

%% Fractions are created by providing the numerator and denominator
FRAC = fr
%%
FRAC = fr(1)
%%
FRAC = fr(1,3)
%%
FRAC = fr(4,3)
%%
FRAC = fr(1,1,3)
%%
FRAC = fr(-4,3)

%% Non-integers are automatically converted (including NaN and Inf)
% keep in mind the default 1e-6 relative accuracy of rat
FRAC = fr(1/3)
%%
FRAC = fr(pi)
%%
FRAC = fr(eps)

%% Fractions can be infinite
FRAC = fr(NaN)
%%
FRAC = fr(0,0)
%%
FRAC = fr(1,0)
%%
FRAC = fr(Inf)
%%
FRAC = fr(-1,0)
%%
FRAC = fr(-Inf)

%% For very large numerators or denominators, use vpi numbers
% doubles can only represent integers exactly up to 2^53-1.
% If you have John D'Errico's Variable Precision Toolbox, use vpi numbers.
% Although fr takes care to avoid integer overflow whenever possible, it is
% recommended that you use vpi numbers.
FRAC = fr(1,2^100)
%%
try
  FRAC = fr(1,vpi(2)^100)
catch
  % need VPI toolbox
end;

%% Arithmetic manipulation of fractions
FRAC = fr(1,3)
%%
FRAC + 1
%%
FRAC + fr(2,3)
%%
FRAC*3
%%
FRAC/4
%%
FRAC^2

%% Work with arrays of fractions
A = fr(1,1:10)
%%
B = A.^2

%%
% Compute the sum
sum(A)

%%
% Or the product
prod(B)

%%
% Convolution
conv(A,B)

%% Using vpi fractions for larger arrays is recommended
%sum(fr(1,1:100))
try
  sum(fr(vpi(1),1:100))
catch
  % need VPI toolbox
end;


%% Relational operators
% All of the standard operators are provided, <, >, <=, >=, ==, ~=
a=fr(1,2)
%%
b=fr(2,3)
%%
1/(a*b) == 3
%% 
a > b
%% 
a <= 1


%% Concatenation using []
[3 2]*[a,b;b,a]*[a^3 b^4]'


%% Conversion from fraction back to double
double(fr(1,3))


%% Extracting the numerator and denominator
[n,d]=rat(fr(4,3))
%%
[k,n,d]=rat(fr(4,3))


%% Extracting digits after the decimal point (including repeating sequence)
% base 10: 1/7=0.142857(142857)
[dig,rep]=digits(fr(1,7))
%%
% base 7: 1/7=0.1(0)
[dig,rep]=digits(fr(1,7),[],7)


%% Extracting a fixed number of digits after the decimal point
% base 10: 1/7=0.142+(3/3500)
[dig,rep]=digits(fr(1,7),3)
%%
% base 7: 1/7=0.100+(0)
[dig,rep]=digits(fr(1,7),3,7)
%%
fr(142,1000)+fr(3,3500)


%% Partial fractions
% partial expansion of 1/12 as sum(a/p^k) where 0<a<p
FPART=partial(fr(1,12))
%%
sum(FPART)


%%
% partial expansion of 1/12 as sum(a/p^k) where 0<a<p^k
FPART=partial(fr(1,12),false)
%%
sum(FPART)


%% Continued fraction expansions and best rational approximations
% continued fraction representation of a fraction:
cf=cfrac(fr(3,8))

%%
% continued fraction expansion of a square root (including repeating term)
[cf,rep]=cfracsqrt(fr(13))


%%
% best rational approximations with a given denominator limit
[r1,r2] = bestrat(fr(cf),rep,100) 


%% Factoring numerator and denominator
F=factor(fr(35,24))
%%
prod(F)


%% Solving linear equations
% use of fractions eliminates many issues due to machine precision
A=fr(ones(2),[2,3;5,7]);
B=fr(ones(2,1),[11;13]);
%%
rref([A,B])
%%
A\B
%%
double(A)\double(B)-double(A\B)

%% Solving rank-deficient linear equations
% fr/mldivide ("\") can handle singular and non-square A
A=[1,1;0,0];
B=[1;0];
%%
% the built-in version does not always find a particular solution
A\B
%%
% the 1-output version of "\" gives NaN values to indicate when the system
% has multiple solutions however substituting 0 for NaN still gives a solution
X=fr(A)\fr(B)
X(isnan(X))=0;
%%
A*X-B
%%
% the 2-output version of "\" provides the general solution X+N*V where V
% is an arbitrary array or fraction array of compatible size
[X,N]=fr(A)\fr(B)
V=fr(5,13);
%%
A*(X+N*V)-B
%%
% when the equations have no solution, "\" returns an empty matrix and a
% warning (this differs from the built-in version)
% e.g. {x=0,x=1}
fr([1;1])\[0;1]
%%
% to find the least squares solution, use lsq:
lsq(fr([1;1]),[0;1])
%%
% See the help files for more details
