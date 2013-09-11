function [fens,gcells,ys] = axialrotorfinmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, nmodules, layers, modulexs, maglabel, varargin)
% creates a structured FAESOR mesh of one fin of a rotor of an axial-flux
% type machine
%
% Syntax
% 
%[fens,gcells] = axialrotorfinmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, layers, circumpoints)
%
%
% Input
%
% Rs - inner radius of rotor shaft 
% 
% Rbi - outer radius of rotor shaft, and inner radius of the rotor fin back
%   iron
% 
% Rmi - inner radius of magnet positions
% 
% Rmo - outer radius of magnet positions
% 
% Rbo - outer rotor fin radius
% 
% tbi - thickness of the fin
% 
% layers - number of layers the mesh will be made up of (before refinement)
%   in the same dimension as tbi, the rotor fin thickness.
% 
% circumpoints - number of circumferential points the mesh will be made up
%   of (before refinement)
%
%
% Output
%
% fens - FAESOR fenodeset object containing the nodes of the rotor fin mesh
%
% gcells - FAESOR gcellset_H20 oject containing the H20 type cells making
%   up the rotor fin mesh
%
%

    Inputs.BackIronRadialPointsPerM = 15;
    Inputs.MagnetRadialPointsPerM = 15;
    Inputs.zshift = 0;
    Inputs.Subset = nmodules;
    if nmodules == 1
        Inputs.ModuleGap = 0;
    else
        Inputs.ModuleGap = 2*modulexs(1);
    end
    
    Inputs = parseoptions(Inputs, varargin);
    
    if nargin < 9
        maglabel = 0;
    end
    
    modulegap = Inputs.ModuleGap;
    
    mergetol = Rs * min([modulegap/10, diff(modulexs)./2, tbi/10]);
    
    % to make rotor mesh, 

    % first make modules
    % 1. make set of Q4 quads of approriate spacing
    % construct spacings
    xs = linspace(0, 2*pi, nmodules+1);
    
    ys = unique( [ Rs, ...
                   linspace(Rbi, Rmi, max(2, round((Rmi-Rbi)*Inputs.BackIronRadialPointsPerM))), ...
                   linspace(Rmi, Rmo, max(2, round((Rmo-Rmi)*Inputs.MagnetRadialPointsPerM))), ...
                   linspace(Rmo, Rbo, max(2, round((Rbo-Rmo)*Inputs.BackIronRadialPointsPerM)))
                 ] );
    
%     if nmodules == 1
%         modulexs = [0, modulexs, 2*pi];
%     end
    
    % create a set of 2D cells we can extrude into the 3D mesh
    [fens, gcells] = Q4_blockx(modulexs, ys, []); 
    
    thismodulexs = modulexs + xs(1);
    
    for i = 2:Inputs.Subset
        
        thismodulexs = modulexs + xs(i);
        
        % create a set of 2D cells we can extrude into the 3D mesh
        [newfens, newgcells] = Q4_blockx(thismodulexs, ys, []);
        
        [fens,gcells,newgcells] = merge_meshes(fens, gcells, newfens, newgcells, 2*eps);
        
        gcells = cat(gcells, newgcells);
        
        [linkfens,linkcells] = Q4_blockx([xs(i)-modulegap/2, xs(i)+modulegap/2], [Rs, Rbi], []);
        
        [fens,gcells,linkgcells] = merge_meshes(fens, gcells, linkfens, linkcells, 2*eps);
        
        gcells = cat(gcells, linkgcells);
        
%         % get the nodes at the edges of the rotor shaft
%         conn = [ fenode_select(fens, struct ('nearestto',[xs(i)-modulegap/2, Rs])), ...
%                  fenode_select(fens, struct ('nearestto',[xs(i)+modulegap/2, Rs])), ...
%                  fenode_select(fens, struct ('nearestto',[xs(i)+modulegap/2, Rbi])), ...
%                  fenode_select(fens, struct ('nearestto',[xs(i)-modulegap/2, Rbi])) ];
%         
%         % make a cell to link the new module and the previous one
%         gcells = cat(gcells, gcellset_Q4(struct( 'conn', conn)));
        
    end
    
    if nmodules > 1
        
        % make a final cell to link the last module to the first when rotated
        [newfens, newgcells] = Q4_blockx([thismodulexs(end), thismodulexs(end)+modulegap], ys(1:2), []);
        
        [fens,gcells,newgcells] = merge_meshes(fens, gcells, newfens, newgcells, []);
        
        gcells = cat(gcells, newgcells);
        
    end
    
    % give the appropriate cells the magnet label
    celllist = gcell_select(fens,gcells,struct ('box', [0, 2*pi, Rmi, Rmo], 'inflate', 2*eps));
    
    % set the cells to have the magnet label
    gcells = setlabels(gcells, maglabel, celllist);
    
    % extrude the mesh into the desired number of layers
    [fens, gcells] = Q4_extrude(fens, gcells, layers, tbi, true );
    % drawmesh({fens, gcells});

    % subdivide mesh into 20 nodes instead of 8 nodes
%     [fens,gcells] = H8_to_H20(fens,gcells);
%     [fens,gcells] = H8_to_H20(fens,gcells);

    % get the locations of the nodes in the mesh
    xyz = get (fens,'xyz');

    % 
    [x,y] = pol2cart(-xyz(:,1),xyz(:,2));
    
    xyz = [x, y, xyz(:,3)];
    
    % rotate them into a cylinder shape
%     for i = 1:count(fens)
% 
%         x = xyz(i,1); 
%         y = xyz(i,2); 
%         z = xyz(i,3);
% 
%         xyz(i,:) = [0, y, z];
% %(0.5*Rs)
%         xyz(i,:) = xyz(i,:) * rotmat([0, 0, -x]); % + [0,0,0];
% 
%     end

    % shift middle of fin to origin
    xyz(:,3) = xyz(:,3) - (tbi/2) + Inputs.zshift;

    % replace the coordinates of the nodes
    fens = set(fens, 'xyz', xyz);

    % Merge together the nodes at the fin edges joint
    [fens, gcells] = merge_nodes(fens, gcells, mergetol);

%     drawmesh({fens, gcells});

end