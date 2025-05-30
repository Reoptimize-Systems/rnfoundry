function y = slmeval(x, slm, evalmode, checkinput)
% slmeval: evaluates a Hermite function, its derivatives, or its inverse
%
% Syntax
%
% y = slmeval(x, slm) 
% y = slmeval(..., evalmode)
% y = slmeval(..., checkinput)
%
% Description
%
% slmeval evaluates a Hermite function, its derivatives, or its inverse, in
% the form produced by slm_engine. See slm_engine, slmset and slm_tutorial
% for more information.
%
% usage 1: y = slmeval(x,slm)           % evaluates y = f(x)
% usage 2: y = slmeval(x,slm,evalmode)  % general form
%
% As opposed to extrapolation as I am, slmeval will not extrapolate. If you
% wanted to extrapolate, then you should have built the spline differently.
% 
% Input
%
%  x - array (any shape) of data to be evaluated through model
%
%    All points are assumed to lie inside the range of the knots. Any which
%    fall outside the first or last knot will be assigned the value of the
%    function at the corresponding knot. NO extrapolation is allowed. (If
%    you really need to extrapolate, then you should have done it when you
%    built the model in the first place. An alternative, if I ever write
%    it, is an extrapolator code, that builds a new, extrapolated
%    function.)
%
%  slm - a shape language model structure, normally constructed by slmfit 
%    or slmengine.
%
%    The fields in this struct are:
%      
%      slm.type : either 'cubic', 'linear', or 'constant'
%
%      slm.knots : a vector of knots (distinct & increasing) There must be
%        at least two knots.
%
%      slm.coef : an array of coefficients for the hermite function.
%
%    If the function is linear or constant Hermite, then slm.coef will have
%    only one column, composed of the value of the function at each
%    corresponding knot. The 'constant' option uses a function value at
%    x(i) to apply for x(i) <= x < x(i+1).
%
%    If the function is cubic Hermite, then slm.coef will have two columns.
%    The first column will be the value of the function at the
%    corresponding knot, the second column will be the correspodign first
%    derivative.
%
% evalmode - (OPTIONAL) numeric flag - specifies what evaluation is to be 
%   done.
%
%   DEFAULT VALUE: 0
% 
%   == 0 --> evaluate the function at each point in x.
%   == 1 --> the first derivative at each point in x.
%   == 2 --> the second derivative at each point in x.
%   == 3 --> the third derivative at each point in x.
%   == -1 --> evaluate the inverse of the function at
%    each point in x, thus y is returned such that x=f(y)
% 
%   Note 1. Piecewise constant functions will return zero for all order
%   derivatives, since I ignore the delta functions at each knot. The
%   inverse operation is also disabled for constant functions.
% 
%   Note 2. Linear hermite functions will return zero for second and higher
%   order derivatives. At a knot point, while technically the derivative is
%   undefined at that location, the slope of the segment to the right of
%   that knot is returned.
% 
%   Note 3. Inverse computations will return the LEFTMOST zero (closest to
%   -inf) in the event that more than one solution exists.
% 
%   Note 4: Inverse of points which fall above the maximum or below the
%   minimum value of the function will be returned as a NaN.
%
% checkinput - true/false flag indicating whether to do input checking.
%   This can be turned off to speed up evaluation. Default is true if not
%   supplied.
%
% Arguments: (output)
%  y       - Evaluated result, the predicted value, i.e., f(x)
%            or f'(x), f''(x), f'''(x), or the functional inverse
%            such that x = f(y). y will have the same shape and
%            size as x.
%
%
% See also: slmpar, ppval, slmengine
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 6/10/09

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
    
    persistent usemexslmeval
    
    if isempty(usemexslmeval)
        if exist ('mexslmeval', 'file') == 3
            usemexslmeval = true;
        else
            usemexslmeval = false;
        end
    end
    
    if nargin < 4
        
        checkinput = true;
        
        % default evalmode is 0, i.e., evaluate f(x)
        if nargin < 3
            evalmode = 0;
        else
            % default evalmode is 0, i.e., evaluate f(x)
            if isempty(evalmode)
                evalmode = 0;
            end
        end
        
    end
    
    if checkinput
        % check for the minimum required information
        if ~isstruct(slm) || ~isfield(slm,'knots') || ...
                ~isfield(slm,'coef') || ~isfield(slm,'degree')
            error('slm must be a struct, with "knots" "coef" and "degree" fields')
        end
    end
    
    % get the shape of x, so we can restore it at the end
    sizex = size(x);
    % then make x a column vector for the internal calculations
    x = x(:);
    % get the knots and make them a column vector
    knots = slm.knots(:);
    % count the number of knots
    nk = length(knots);

    % evaluate the appropriate class of piecewise hermite function
    switch slm.degree
        case 1
            nx = numel(x);
            dx = diff(knots);
            if checkinput
                % verify that coef is the right size
                sc=size(slm.coef);
                if (sc(1) ~= nk)
                    error('Improper size of coef field for these knots')
                elseif (sc(2)~=1)
                    error('Improper size of coef field for a linear Hermite')
                end
            end

            % which knot interval does each point fall in?
            % only treat evalmode == 0 or 1 right now
            if ismember(evalmode,[0 1])
                % clip the points first
                x = min(knots(end),max(x,knots(1)));

                % and use histc to bin the points.
                [~,xbin] = histc(x,knots);

                % any point which falls at the top end, is said to
                % be in the last bin.
                xbin(xbin==nk)=nk-1;

            else
                % we need to treat the inverse carefully, and we
                % can ignore where a point falls for any higher
                % derivatives requested. do this part below.
                xbin = zeros(nx,1);

            end

            % Evaluation mode:
            switch evalmode
                case 0
                    % f(x)
                    t = (x-knots(xbin))./dx(xbin);
                    y = slm.coef(xbin) + (slm.coef(xbin+1)-slm.coef(xbin)).*t;

                case 1
                    % first derivative
                    y = (slm.coef(xbin+1)-slm.coef(xbin))./dx(xbin);

                case {2 3}
                    % higher derivatives of a piecewise linear function
                    % are all zero
                    y = zeros(sizex);

                case -1
                    % functional inverse
                    % here the problem of which bin we fall into is not
                    % trivial, since the function need not be monotone.
                    % we might also have constant segments.
                    % Identify the leftmost interval which contains the
                    % point in question.

                    % first, clip the data in terms of function value
                    miny = min(slm.coef);
                    maxy = max(slm.coef);
                    y = repmat(NaN,nx,1);
                    k = find((x<=maxy) & (x>=miny));

                    % determine which knot interval to look in for
                    % each point.
                    ybin = nmbs(x(k),slm.coef);

                    % rule out the divide by zero cases for intervals
                    % where the function was constant
                    L = (slm.coef(ybin+1) == slm.coef(ybin));
                    if any(L)
                        y(k(L)) = knots(ybin(L));
                        ybin(L)=[];
                        k(L)=[];
                    end

                    % inverse interpolation
                    y(k) = knots(ybin) + (x(k) - slm.coef(ybin)).* ...
                        (knots(ybin+1)-knots(ybin))./(slm.coef(ybin+1)-slm.coef(ybin));

                otherwise
                    % anything else
                    error('Evalmode must be one of [0 1 2 3 -1]')
            end

        case 3
            
                
            nx = numel(x);
            dx = diff(knots);
            if checkinput
                % verify shape of coef
                sc=size(slm.coef);
                if (sc(1) ~= nk)
                    error('Improper size of coef field for these knots')
                elseif (sc(2)~=2)
                    error('Improper size of coef field for a cubic Hermite')
                end
            end

            if evalmode == -1
                % we need to treat the inverse more carefully
                xbin = zeros(nx,1);
            else
                if ~usemexslmeval
                    % clip the points first
                    x = min(knots(end),max(x,knots(1)));

                    % and use histc to bin the points.
                    [~,xbin] = histc(x,knots);

                    % any point which falls at the top end is said to
                    % be in the last bin.
                    xbin(xbin==nk)=nk-1;
                end
            end
            
            % Evaluation mode:
            switch evalmode
                case 0
                    if usemexslmeval
                        y = mexslmeval (x, knots, slm.coef, 0);
                    else
                        % f(x)
                        t = (x-knots(xbin))./dx(xbin);
                        t2 = t.^2;
                        t3 = t.^3;
                        s2 = (1-t).^2;
                        s3 = (1-t).^3;
                        y = (-slm.coef(xbin,2).*(s3-s2) + ...
                            slm.coef(xbin+1,2).*(t3-t2)).*dx(xbin) + ...
                            slm.coef(xbin,1).*(3*s2-2*s3) + ...
                            slm.coef(xbin+1,1).*(3*t2-2*t3);
                    end
                case 1
                    if usemexslmeval
                        y = mexslmeval (x, knots, slm.coef, 1);
                    else
                        % first derivative for the cubic case
                        t = (x-knots(xbin))./dx(xbin);
                        t2 = t.^2;
                        s = 1-t;
                        s2 = (1-t).^2;
                        y = -slm.coef(xbin,2).*(-3*s2+2*s) + ...
                            slm.coef(xbin+1,2).*(3*t2-2*t) + ...
                            (slm.coef(xbin,1).*(-6*s+6*s2) + ...
                            slm.coef(xbin+1,1).*(6*t-6*t2))./dx(xbin);
                    end
                case 2
                    if usemexslmeval
                        y = mexslmeval (x, knots, slm.coef, 2);
                    else
                        % second derivative of a cubic
                        t = (x-knots(xbin))./dx(xbin);
                        s = 1-t;
                        y = (-slm.coef(xbin,2).*(6*s - 2) + ...
                            slm.coef(xbin+1,2).*(6*t - 2))./dx(xbin) + ...
                            (slm.coef(xbin,1).*(6 - 12*s) + ...
                            slm.coef(xbin+1,1).*(6 - 12*t))./(dx(xbin).^2);
                    end
                case 3
                    if usemexslmeval
                        y = mexslmeval (x, knots, slm.coef, 3);
                    else
                        % third derivative
                        y = 6*(slm.coef(xbin,2) + slm.coef(xbin+1,2))./(dx(xbin).^2) + ...
                            12*(slm.coef(xbin,1) - slm.coef(xbin+1,1))./(dx(xbin).^3);
                    end
                case -1
                    % functional inverse
                    % here the problem of which bin we fall into is not
                    % trivial, since the function need not be monotone.
                    % we might also have constant segments.

                    % first, convert the spline into a pp form.
                    pp = slm2pp(slm);
                    coefs = pp.coefs;

                    % scale the cubic polys so they live on [0,1].
                    % this will stabilize things a bit, as well as make
                    % the tests easier later on. I could do this in one
                    % line with a variety of tricks, but feeling lazy...
                    coefs(:,1) = coefs(:,1).*(dx.^3);
                    coefs(:,2) = coefs(:,2).*(dx.^2);
                    coefs(:,3) = coefs(:,3).*dx;

                    % Identify the leftmost interval which contains the
                    % point in question. The problem is the cubic segments
                    % might not be monotone. nmbs will try though.
                    binedges = slm.coef(:,1);
                    ybin = nmbs(x,binedges);

                    % a simple solution is to test every interval up to
                    % and including the interval found in ybin, but no
                    % further. We can stop there since we know a solution
                    % must exist in that interval. If no interval was found,
                    % then ybin will be a NaN. In that event, we must search
                    % every interval, since the curve need not be monotone.
                    ybin(isnan(ybin) | (ybin == nk)) = nk - 1;

                    % just a loop over roots here. not terribly efficient,
                    % but I don't terribly want to code a vectorized cubic
                    % solver for this problem.
                    y = NaN(size(x));
                    tol = 1000*eps;
                    for j = 1:nx
                        flag = true;
                        i = 1;
                        while flag && (i <= ybin(j))
                            Ci = coefs(i,:);

                            % offset the constant term. a root of this
                            % cubic is what we want.
                            Ci(4) = Ci(4) - x(j);

                            % scale the poly coefficients to improve things
                            % yet some more.
                            Ci = Ci./max(abs(Ci(1:3)));

                            % get the roots. at least one must be real.
                            Ri = roots(Ci);

                            % of the real roots, is one of them strictly
                            % inside [0,1]? If more than one is, take the
                            % smallest root.
                            k = (abs(imag(Ri)) == 0) & (real(Ri) >= 0) & (real(Ri) <= 1);
                            if any(k)
                                k = find(k,1,'first');
                                y(j) = knots(i) + dx(i)*real(Ri(k));
                                flag = false;
                                continue
                            end

                            % if we did not succeed in the last test,
                            % of the real roots, is one of them within
                            % inside [-tol,1+tol]? If more than one is,
                            % take the smallest root.
                            k = (abs(imag(Ri)) == 0) & (real(Ri) >= -tol) & (real(Ri) <= (1 + tol));
                            if any(k)
                                k = find(k,1,'first');
                                y(j) = knots(i) + dx(i)*real(Ri(k));
                                flag = false;
                                continue
                            end

                            % increment i, check the next interval.
                            i = i + 1;
                        end
                    end

                otherwise
                    % anything else
                    error('Evalmode must be one of [0 1 2 3 -1]')
            end

        case 0
            % piecewise constant function (discontinuous at the knots)
            % verify that coef is the right size
            if checkinput
                sc=size(slm.coef);
                if (sc(1) ~= (nk-1))
                    error('Improper size of coef field for these knots')
                elseif (sc(2)~=1)
                    error('Improper size of coef field for a piecewise constant function')
                end
            end

            % Evaluation mode:
            switch evalmode
                case 0
                    % f(x)
                    if numel(knots) > 2
                        % clip the points first
                        x = min(knots(end),max(x,knots(1)));

                        % which knot interval does each point fall in?
                        % and use histc to bin the points.
                        [~,xbin] = histc(x,knots);

                        % any point which falls at the top end, is said to
                        % be in the last bin.
                        xbin(xbin==nk)=nk-1;
                        
                        y = slm.coef(xbin);
                    else
                        % only two knots, so we are in bin 1 at all points
                        y = repmat(slm.coef(1), size(x));
                    end

                case {1 2 3}
                    % all derivatives of a piecewise constant function
                    % will be zero if we ignore the delta functions
                    % at each knot.
                    y = zeros(sizex);

                case -1
                    % functional inverse - no inverse exists for a
                    % discontinuous function.
                    y = reshape(NaN,sizex);

                otherwise
                    % anything else
                    error('Evalmode must be one of [0 1 2 3 -1]')
            end

        otherwise
            error('slmeval only handles the constant, linear, or cubic Hermite case')
    end

    % when all done, make sure the size of the output is the
    % same as the input was.
    y = reshape(y,sizex);

end


% ===================================================
% ================ begin subfunction ================
% ===================================================
function bind = nmbs(x,binedges)
% nmbs: non-monotone bin search
%
% Finds the leftmost bin that contains each x. Assumes
% that x has already been clipped so it must lie in
% some bin.
nb = length(binedges);
nx = length(x);

% if the bins are actually monotone, then just use histc
db = diff(binedges);
if all(db>0)
    % increasing bins
    [~,bind] = histc(x,binedges);
    bind(bind==nb)=nb-1;

elseif all(db<0)
    % decreasing sequence of edges
    [~,bind] = histc(-x,-binedges);
    bind(bind==nb)=nb-1;

else
    % non-monotone sequence of edges. Do this one the
    % hard way. Find the first bin that fits.
    i = 1;
    j = 1:nx;
    bind = ones(nx,1);
    while ~isempty(j) && (i<nb)
        k = ((binedges(i)>=x(j)) & (binedges(i+1)<=x(j))) | ...
            ((binedges(i)<=x(j)) & (binedges(i+1)>=x(j)));

        if any(k)
            bind(j(k)) = i;
            j(k)=[];
        end

        % increment bin we will look in
        i=i+1;

    end

end

end






