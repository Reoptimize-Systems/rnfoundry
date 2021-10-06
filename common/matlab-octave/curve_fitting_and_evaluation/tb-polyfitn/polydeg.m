function n = polydeg(x,y)
%POLYDEG Find an optimal degree for polynomial fitting.
%   N = POLYDEG(X,Y) finds the optimal degree for polynomial fitting
%   according to the Akaike's information criterion.
%
%   Assuming that you want to find the degree N of a polynomial that fits
%   the data Y(X) best in a least-squares sense, the Akaike's information
%   criterion is defined by:
%       2*(N+1) + n*(log(2*pi*RSS/n)+1)
%   where n is the number of points and RSS is the residual sum of squares.
%   The optimal degree N is defined as that which minimizes <a
%   href="matlab:web('http://en.wikipedia.org/wiki/Akaike_Information_Criterion')">AIC</a>.
%
%   Notes:
%   -----
%   If the number of data is small, POLYDEG may tend to return:
%   N = (number of points)-1.
%
%   ORTHOFIT (FEX #16574) is more appropriate than POLYFIT for polynomial
%   fitting with relatively high degrees (available <a
%   href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=16574')">here</a>).
%
%   Examples:
%   --------
%   load census
%   n = polydeg(cdate,pop)
%
%   x = linspace(0,10,300);
%   y = sin(x.^3/100).^2 + 0.05*randn(size(x));
%   n = polydeg(x,y)
%   ys = orthofit(x,y,n);
%   plot(x,y,'.',x,ys,'k')
%
%   Damien Garcia, 02/2008, revised 07/2008
%
%   See also POLYFIT, ORTHOFIT.

%% Check input arguments
% ---
%
error(nargchk(2,2,nargin));

x = x(:);
y = y(:);
N = length(x);

if ~isequal(N,length(y))
    error('MATLAB:polydeg:XYNumelMismatch',...
          'X and Y must have same number of elements.')
elseif ~isfloat(x) || ~isfloat(y)
    error('MATLAB:polydeg:NonFloatingXY',...
    'X and Y must be floating point arrays.')
end

%% Search the optimal degree minimizing the Akaike's information criterion
% ---
%  y(x) are fitted in a least-squares sense using a polynomial of degree n
%  developed in a series of orthogonal polynomials.
%

% Turn warning messages off
warn01 = warning('query','MATLAB:log:logOfZero');
warn02 = warning('query','MATLAB:divideByZero');
warning('off','MATLAB:log:logOfZero')
warning('off','MATLAB:divideByZero')

% 0th order
p = mean(y);
ys = ones(N,1)*p;
AIC = 2+N*(log(2*pi*sum((ys-y).^2)/N)+1)+...
    4/(N-2);

p = zeros(2,2);
p(1,2) = mean(x);
PL = ones(N,2);
PL(:,2) = x-p(1,2);

n = 1;
nit = 0;

% While-loop is stopped when a minimum is detected. 3 more steps are
% required to take AIC noise into account and to ensure that this minimum
% is a non-local minimum.
% ---
while nit<3
    
    % -- Orthogonal polynomial fitting (see also ORTHOFIT)
    if n>1
        p(1,n+1) = sum(x.*PL(:,n).^2)/sum(PL(:,n).^2);
        p(2,n+1) = sum(x.*PL(:,n-1).*PL(:,n))/sum(PL(:,n-1).^2);
        PL(:,n+1) = (x-p(1,n+1)).*PL(:,n)-p(2,n+1)*PL(:,n-1);
    end
    
    tmp = sum(repmat(y,[1,n+1]).*PL)./sum(PL.^2);
    ys = sum(PL.*repmat(tmp,[N,1]),2);
    
    % -- Akaike's Information Criterion
    aic = 2*(n+1)+N*(log(2*pi*sum((ys-y(:)).^2)/N)+1)+...
        2*(n+1)*(n+2)/(N-n-2); % correction for small sample sizes
    
    if aic>=AIC
        nit = nit+1;
    else
        nit = 0;
        AIC = aic;
    end
    n = n+1;
        
    if n>=N, break, end
    
end

n = n-nit-1;

% Go back to previous warning states
warning(warn01)
warning(warn02)

