function [fens,gcells] = axialrotorshaftmesh3d(Rs, Rbi, tshaftsec, nmodules, layers, modulexs, varargin)

    Inputs.Subset = nmodules;
    if nmodules == 1
        Inputs.ModuleGap = 0;
    else
        Inputs.ModuleGap = 2*pi/(100*nmodules);
    end
    
    Inputs = parseoptions(Inputs, varargin);
    
    % first make modules
    % 1. make set of Q4 quads of approriate spacing
    % construct spacings
    xs = linspace(0, 2*pi, nmodules+1);
        
    modulegap = Inputs.ModuleGap;
    
    ys = [ Rs, Rbi ];
    
    % create a set of 2D cells we can extrude into the 3D mesh
    [fens, gcells] = Q4_blockx(modulexs, ys, []); 
    
    thismodulexs = modulexs + xs(1);
    
    for i = 2:Inputs.Subset
        
        thismodulexs = modulexs + xs(i);
        
        % create a set of 2D cells we can extrude into the 3D mesh
        [newfens, newgcells] = Q4_blockx(thismodulexs, ys, []);
        
        [fens,gcells,newgcells] = merge_meshes(fens, gcells, newfens, newgcells, 2*eps);
        
        gcells = cat(gcells, newgcells);
        
        [linkfens,linkgcells] = Q4_blockx([xs(i)-modulegap/2, xs(i)+modulegap/2], ys, []);
        
        [fens,gcells,linkgcells] = merge_meshes(fens, gcells, linkfens, linkgcells, 2*eps);
        
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
        [newfens, newgcells] = Q4_blockx([thismodulexs(end), thismodulexs(end)+modulegap], ys, []);
        
        [fens,gcells,newgcells] = merge_meshes(fens, gcells, newfens, newgcells, 2*eps);
        
        gcells = cat(gcells, newgcells);
        
    end
    
    % extrude the mesh into the desired number of layers
    [fens, gcells] = Q4_extrude(fens, gcells, layers, @(xyz,k)[xyz(1:2),k*(tshaftsec/layers)], true );
    % drawmesh({fens, gcells});

    % subdivide mesh into 20 nodes instead of 8 nodes
%     [fens,gcells] = H8_to_H20(fens,gcells);

    % get the locations of the nodes in the mesh
    xyz = get (fens,'xyz');
    
    [x,y] = pol2cart(-xyz(:,1),xyz(:,2));
    
    xyz = [x, y, xyz(:,3)];

    % shift middle of fin to origin
    xyz(:,3) = xyz(:,3) - (tshaftsec/2);

    % replace the coordinates of the nodes
    fens = set(fens, 'xyz', xyz);

    % Merge together the nodes at the fin edges joint
    [fens, gcells] = merge_nodes(fens, gcells, 2*eps);

end