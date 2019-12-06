
R = 0.25;
n = 25;
cut_z_max = 0.25 * R;

[vertices, faces, info, hfig, hax] = mesh.quad.hemisphere (R, n, 'CutZMax', cut_z_max, 'Draw', false);

figure;
patch('vertices', vertices', 'faces', faces', 'FaceAlpha', 0, 'EdgeColor', 'k');

hold on
scatter3 (vertices(1,info.CutTopInds), vertices(2,info.CutTopInds), vertices(3,info.CutTopInds), 'xb');
scatter3 (vertices(1,info.EndInds), vertices(2,info.EndInds), vertices(3,info.EndInds), 'or');
hold off
        
axis equal




%%
% 
% R = 0.25;
% n = 25;
% cut_z_max = 1.1 * R;
% 
% [vertices, faces, info, hfig, hax] = mesh.quad.hemisphere (R, n, 'CutZMax', cut_z_max, 'Draw', true);
% 
% figure;
% patch('vertices', vertices', 'faces', faces', 'FaceAlpha', 0, 'EdgeColor', 'k');
% 
% hold on
% scatter3 (vertices(1,info.CutTopInds), vertices(2,info.CutTopInds), vertices(3,info.CutTopInds), 'xb');
% scatter3 (vertices(1,info.EndInds), vertices(2,info.EndInds), vertices(3,info.EndInds), 'or');
% hold off
%         
% axis equal