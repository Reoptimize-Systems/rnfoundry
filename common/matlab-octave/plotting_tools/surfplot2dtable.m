function [newx, newy] = surfplot2dtable (x, y, data, varargin)
% makes a surf plot of data in the form of a table
%
% Syntax
%
% [newx, newy] = surfplot2dtable (x, y, data)
% [...] = surfplot2dtable (..., 'Parameter', value)
%
% Description
%
% surfplot2dtable makes a surf plot of data provided in the form of a 2d
% matrix, with associated values for the rows (x values) and columns (y
% values).
%
% Input
%
%  x - vector of m values the same length as the number of rows in 'data'.
%
%  y - vector of n values the same length as the number of columns in
%    'data'.
%
%  data - (m x n) matrix of numeric values. A 3D triangulated surface plot
%    will be created based on this data and the corresponding x and y
%    values.
%
% Output
%
%  newx - x coordinates at which each point in 'data' is plotted
%
%  newy - y coordinates at which each point in 'data' is plotted
%
%
% See Also:  contour2dtable.m contourf2dtable.m
%


    options.XLabel = '';
    options.YLabel = '';
    options.ZLabel = '';
    options.Title = '';
    options.TriMeshArgs = {};
    options.Axes = [];
    
    options = parse_pv_pairs (options, varargin);

    [newx, newy, newz] = makevars2dtable23dplot(x, y, data);
    
    clear data;
    
    % create a Delaunay triangulation.
    tri = delaunay(newx,newy);

    if isempty (options.Axes)
        figure;
        axes;
    end
    
    % Plot it with TRISURF
    trimesh(tri, newx, newy, newz, options.TriMeshArgs{:});
    axis vis3d

    % Clean it up

%     axis off
%     l = light('Position',[-50 -15 29])
%     set(gca,'CameraPosition',[208 -50 7687])
%     lighting phong
%     shading interp
%     colorbar EastOutside

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