%%
n = 16;
R = 0.25;
[Vertices, Faces, edge_inds] = mesh.quad.spherepatchgrid (R, n, ...
                                                'ProjectionType', 'equidistance', ...
                                                'CutZMax', 0.2 * R );


fig = figure(1); clf(fig);
pos=get(fig, 'Position');
pos([3 4]) = 400;
set(fig, 'Position', pos);
ax=axes('Parent',fig,'pos',[0 0 1 1]);
patch('Vertices', Vertices', 'Faces', Faces', 'FaceAlpha', 0, 'parent',ax);
axis(ax,'equal')
% axis(ax,'off')
view(ax,3)

hold on
scatter3 (Vertices(1, edge_inds), Vertices(2, edge_inds), Vertices(3, edge_inds), 'or');
hold off



