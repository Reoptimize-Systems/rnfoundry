function [xq,yq,zq,C,hc,hax,hfig] = contour2dtable(x, y, data, varargin)
% creates a contour plot of the data in a 2D table
%
% Syntax
%
% contour2dtable(x, y, data, ...)
%
% Inputs
% 
%  x - vector of values representing the independent variable variables of
%   the rows of the table, must be linearly spaced
%
%  y - vector of values representing the independent variable variables of
%   the columns of the table, must be linearly spaced
%
%  data - (2 x 2) matrix of data containing the dependent variables
%    corresponding to each combination of the independent values provided
%    in the x and y vectors
%
% Examples
%
% x = 1:10
% y = 0.5:0.5:5
% data = bsxfun(@times, x, y')
% contour2dtable(x, y, data)
%
% See also: contourf2dtable.m surfplot2dtable.m
%

% Copyright Richard Crozier 2015

    options.XLabel = '';
    options.YLabel = '';
    options.ZLabel = '';
    options.Title = '';
    options.ContourArgs = {};
    options.Axes = [];
    
    options = parse_pv_pairs (options, varargin);
    
    [newx, newy, newz] = makevars2dtable23dplot (x, y, data);
    
    % now mesh the data
    [xq,yq] = meshgrid (x,y);
    
    zq = griddata (newx, newy, newz, xq, yq);
    
    if isempty (options.Axes)
        hfig = figure;
        hax = axes;
    else
        hax = options.Axes;
        hfig = get (hax, 'Parent');
    end
    
    % and create the mesh plot
    [C,hc] = contour (hax, xq, yq, zq, options.ContourArgs{:});

    if ~isempty (options.XLabel)
        xlabel (options.XLabel);
    end

    if ~isempty (options.XLabel)
        ylabel (options.YLabel);
    end
    
    if ~isempty (options.XLabel)
        zlabel (options.ZLabel);
    end
    
    if ~isempty (options.Title)
        title (options.Title);
    end
    
end