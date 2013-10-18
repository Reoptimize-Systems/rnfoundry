function [newx, newy] = surfplot2dtable(x, y, data, varargin)

    [newx, newy, newz] = makevars2dtable23dplot(x, y, data);
    
    clear data;
    
    % create a Delaunay triangulation.
    tri = delaunay(newx,newy);

    % Plot it with TRISURF
    trimesh(tri, newx, newy, newz, varargin{:});
    axis vis3d

    % Clean it up

%     axis off
%     l = light('Position',[-50 -15 29])
%     set(gca,'CameraPosition',[208 -50 7687])
%     lighting phong
%     shading interp
%     colorbar EastOutside

end