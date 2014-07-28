function [mhs, yy] = monocubicherm (varargin)
% fits a monotonic cubic hermite spline btween two point with the specified
% slopes at the start and end
%
% Syntax
% 
% [mhs] = monocubicherm (x1, y1, x2, y2, dy1, dy2)
% [mhs, yy] = monocubicherm (x1, y1, x2, y2, dy1, dy2, xx)
% [mhs, yy] = monocubicherm (mhs, xx)
%
% 

    % create the fit

    if nargin >= 6
        
        mhs.x1 = varargin{1};
        mhs.y1 = varargin{2};    
        mhs.x2 = varargin{3};
        mhs.y2 = varargin{4};
        mhs.dy1 = varargin{5};
        mhs.dy2 = varargin{6};

        mhs.si = (mhs.y2 - mhs.y1) / (mhs.x2 - mhs.x1);

        mhs.ci = (3 .* mhs.si - 2 .* mhs.dy1 - mhs.dy2) ./ (mhs.x2 - mhs.x1);

        mhs.di = (mhs.dy1 + mhs.dy2 - 2 .* mhs.si) / (mhs.x2 - mhs.x1).^2;
    
        if nargin == 7
            xx = varargin {7};
        else
            xx = [];
        end
        
    elseif nargin == 2
        mhs = varargin{1};
        xx = varargin{2};
    else
        error ('invalid number of arguments.')
    end

    if ~isempty (xx)
        yy = mhs.y1 + mhs.dy1 .* (xx - mhs.x1) + mhs.ci .* (xx - mhs.x1).^2 + mhs.di .* (xx - mhs.x1).^3;
    else
        yy = [];
    end
    
end

