function [fens,gcells] = axialinnerrotorshaftmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, outersep, nmodules, nstages, layers, circumpoints)

    
    % first make modules
    % 1. make set of Q4 quads of approriate spacing
    % construct spacings
    xs = linspace(0, 2*pi, nmodules+1);
    
    if nmodules == 1
        modulegap = 0;
    else
        modulegap = 2*pi/(100*nmodules);
    end
    
    modulexs = linspace(modulegap/2, 2*pi/nmodules - modulegap/2, ceil(circumpoints/nmodules));
    
    bipointsperm = 10;
    
    ys = unique( [ Rs, ...
                   linspace(Rbi, Rmi, max(2, round((Rmi-Rbi)*bipointsperm))), ...
                 ] );
    
    % create a set of 2D cells we can extrude into the 3D mesh
    [fens, gcells] = Q4_blockx(modulexs, ys, []); 
    
    for i = 2:nmodules
        
        thismodulexs = modulexs + xs(i);
        
        % create a set of 2D cells we can extrude into the 3D mesh
        [newfens, newgcells] = Q4_blockx(thismodulexs, ys, []);
        
        [fens,gcells,newgcells] = merge_meshes(fens, gcells, newfens, newgcells, 0);
        
        gcells = cat(gcells, newgcells);
        
        % get the nodes at the edges of the rotor shaft
        conn = [ fenode_select(fens, struct ('nearestto',[xs(i)-modulegap/2, Rs])), ...
                 fenode_select(fens, struct ('nearestto',[xs(i)+modulegap/2, Rs])), ...
                 fenode_select(fens, struct ('nearestto',[xs(i)+modulegap/2, Rbi])), ...
                 fenode_select(fens, struct ('nearestto',[xs(i)-modulegap/2, Rbi])) ];
        
        % make a cell to link the new module and the previous one
        gcells = cat(gcells, gcellset_Q4(struct( 'conn', conn)));
        
    end
    
    % give the appropriate cells the magnet label
    celllist = gcell_select(fens,gcells,struct ('box', [0, 2*pi, Rmi, Rmo]));
    
    % set the cells to have the magnet label
    gcells = setlabels(gcells, maglabel, celllist);
    
    % extrude the mesh into the desired number of layers
    [fens, gcells] = Q4_extrude(fens, gcells, layers, @(xyz,k)[xyz(1:2),k*(tbi/layers)], true );
    % drawmesh({fens, gcells});

    % subdivide mesh into 20 nodes instead of 8 nodes
    [fens,gcells] = H8_to_H20(fens,gcells);

    % get the locations of the nodes in the mesh
    xyz = get (fens,'xyz');

    % rotate them into a cylinder shape
    for i = 1:count(fens)

        x = xyz(i,1); 
        y = xyz(i,2); 
        z = xyz(i,3);

        xyz(i,:) = [0, y-(0.5*Rs), z];

        xyz(i,:) = xyz(i,:) * rotmat([0, 0, -x]) + [0,(0.5*Rs),0];

    end

    % shift middle of fin to origin
    xyz(:,3) = xyz(:,3) - (tbi/2);

    % replace the coordinates of the nodes
    fens = set(fens, 'xyz', xyz);

    % Merge together the nodes at the fin edges joint
    [fens, gcells] = merge_nodes(fens, gcells, 2*eps);

end