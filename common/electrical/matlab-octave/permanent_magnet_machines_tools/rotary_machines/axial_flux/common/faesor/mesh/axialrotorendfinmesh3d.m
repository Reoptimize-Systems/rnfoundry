function [fens,gcells,modulexs,modulegap] = axialrotorendfinmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, tsuppb, tausupp, nmodules, nmodulesupports, whichend, disclayers, supportlayers, circumpoints, maglabel, varargin)
% creates a structured FAESOR mesh of one fin of a rotor of an axial-flux
% type machine
%
% Syntax
% 
% [fens,gcells] = axialrotorendfinmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, tsuppb, tausupp, nmodules, nmodulesupports, disclayers, supportlayers, circumpoints)
% [fens,gcells] = axialrotorendfinmesh3d(..., maglabel)
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
% disclayers - number of disclayers the mesh will be made up of (before refinement)
%   in the same dimension as tbi, the rotor fin thickness.
% 
% circumpoints - number of circumferential points the mesh will be made up
%   of (before refinement)
%
% maglabel - (optional) label for cells in the magnet region
%
% Output
%
% fens - FAESOR fenodeset object containing the nodes of the rotor fin mesh
%
% gcells - FAESOR gcellset_H20 oject containing the H20 type cells making
%   up the rotor fin mesh
%
%

    options.BackIronRadialPointsPerM = 15;
    options.MagnetRadialPointsPerM = 15;
    options.Subset = nmodules;
    if nmodules == 1
        options.ModuleGap = 0;
    else
        options.ModuleGap = 2*pi/(100*nmodules);
    end
    
    options = parseoptions(options, varargin);
    
    if nargin < 13
        maglabel = 0;
    end
    
    % to make rotor mesh,
    
    % first make modules
    % 1. make set of Q4 quads of approriate spacing
    % construct spacings
    xs = linspace(0, 2*pi, nmodules+1);
    
    % convert rib width to an angular width
    tausupp = tausupp / Rbi;
    
    modulegap = options.ModuleGap;

    if strncmpi(whichend, 'top', 1)
        isbottom = false;
    elseif strncmpi(whichend, 'bottom', 1)
        isbottom = true;
    else
        error('whichend must be ''top'' or ''bottom'', or some unambiguous substring')
    end
    
    if nmodulesupports >= 1
        
        suppcenters = linspace(modulegap/2, 2*pi/nmodules - modulegap/2, nmodulesupports+2);
        
        suppcenters = suppcenters(2:end-1);
        
        suppxs = [suppcenters-tausupp/2, suppcenters+tausupp/2];
        
    else
        suppcenters = [];
        suppxs = [];
    end
    
    modulexs = sort(unique([suppxs, linspace(modulegap/2, 2*pi/nmodules - modulegap/2, ceil(circumpoints/nmodules))]));
    
%     mergetol = max(100*eps, min([modulegap/2, diff(modulexs)./2]));
    
    mergetol = Rs * min([modulegap/10, diff(modulexs)./2, tbi/10]);
    
    %     options.zshift = tbi/2;
    
    [finfens,fingcells,ys] = axialrotorfinmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, nmodules, disclayers, modulexs, maglabel, options);
    
    thismodulexs = modulexs + xs(1);
    
    % if there's no supports on the module, we are done, as the disc has
    % been drawn, and there's nothing else to add, otherwise the supports
    % and extra shaft must be drawn and merged with the disc
    if nmodulesupports >= 1 && ~isempty(tsuppb) && (tsuppb > 2*eps(tsuppb)) 
        
        % create cells to make the shaft end part
        [shaftfens, shaftgcells] = Q4_blockx(modulexs, ys(1:2), []);
        
        % create cells to make the end supports
        thissuppxs = modulexs(modulexs >= (suppcenters(1)-tausupp/2) & modulexs <= (suppcenters(1)+tausupp/2));
        
        % create a set of 2D cells we can extrude into the 3D mesh
        [suppfens, suppgcells] = Q4_blockx(thissuppxs, ys(2:end), []);
        
        for i = 2:numel(suppcenters)
            
            thissuppxs = modulexs(modulexs >= (suppcenters(i)-tausupp/2) & modulexs <= (suppcenters(i)+tausupp/2));
            
            % create a set of 2D cells we can extrude into the 3D mesh
            [newsuppfens, newsuppgcells] = Q4_blockx(thissuppxs, ys(2:end), []);
            
            [suppfens,suppgcells,newsuppgcells] = merge_meshes(suppfens, suppgcells, newsuppfens, newsuppgcells, 2*eps);
            
            suppgcells = cat(suppgcells, newsuppgcells);
            
        end
        
        for i = 2:options.Subset
            
            thismodulexs = modulexs + xs(i);
            
            % create cells to make the shaft end part
            [newshaftfens, newshaftgcells] = Q4_blockx(thismodulexs, ys(1:2), []);
            
            [shaftfens,shaftgcells,newshaftgcells] = merge_meshes(shaftfens, shaftgcells, newshaftfens, newshaftgcells, 2*eps);
            
            shaftgcells = cat(shaftgcells, newshaftgcells);

            [linkfens,linkgcells] = Q4_blockx([xs(i)-modulegap/2, xs(i)+modulegap/2], ys(1:2), []);
            
            [shaftfens,shaftgcells,linkgcells] = merge_meshes(shaftfens, shaftgcells, linkfens, linkgcells, mergetol);
            
            shaftgcells = cat(shaftgcells, linkgcells);
            
%             % get the nodes at the edges of the rotor shaft
%             conn = [ fenode_select(shaftfens, struct ('nearestto',[xs(i)-modulegap/2, Rs])), ...
%                 fenode_select(shaftfens, struct ('nearestto',[xs(i)+modulegap/2, Rs])), ...
%                 fenode_select(shaftfens, struct ('nearestto',[xs(i)+modulegap/2, Rbi])), ...
%                 fenode_select(shaftfens, struct ('nearestto',[xs(i)-modulegap/2, Rbi])) ];
%             
%             % make a cell to link the new module and the previous one
%             shaftgcells = cat(shaftgcells, gcellset_Q4(struct( 'conn', conn)));
            
            for ii = 1:numel(suppcenters)
                
                thissuppxs = thismodulexs(modulexs >= (suppcenters(ii)-tausupp/2) & modulexs <= (suppcenters(ii)+tausupp/2));
                
                % create a set of 2D cells we can extrude into the 3D mesh
                [newsuppfens, newsuppgcells] = Q4_blockx(thissuppxs, ys(2:end), []);
                
                [suppfens,suppgcells,newsuppgcells] = merge_meshes(suppfens, suppgcells, newsuppfens, newsuppgcells, mergetol);
                
                suppgcells = cat(suppgcells, newsuppgcells);
                
            end
            
        end
        
        if nmodules > 1
            
            % make final cell to link shaft cylinder
            [linkfens,linkgcells] = Q4_blockx([thismodulexs(end), thismodulexs(end)+modulegap], ys(1:2), []);
            
            [shaftfens,shaftgcells,linkgcells] = merge_meshes(shaftfens, shaftgcells, linkfens, linkgcells, mergetol);
        
            shaftgcells = cat(shaftgcells, linkgcells);
            
        end
        
        
        [shaftfens, shaftgcells] = Q4_extrude(shaftfens, shaftgcells, supportlayers, tsuppb, true );
        
        if isbottom
            shaftfens = translate(shaftfens, [0,0,-(tsuppb+tbi/2)]);
        else
            shaftfens = translate(shaftfens, [0,0,tbi/2]);
        end
            
        if nmodulesupports >= 1 || (tsuppb < 2*eps(tsuppb)) || isempty(tsuppb)
            
            % extrude the mesh to make the shaft
            if isbottom
                extrudedir = -1;
                % reverse the connectivity
                suppconn = getconn(suppgcells);
                suppgcells = setconn(suppgcells, fliplr(suppconn));
            else
                extrudedir = 1;
            end
        
            % create the supports if there are any
            minsuppheight = max(1.1*supportlayers*mergetol, 0.01*tsuppb);
            
            [suppfens, suppgcells] = Q4_extrude(suppfens, suppgcells, supportlayers, ...
                                                @(xyz,k)[xyz(1:2),k*(extrudedir*tsuppbextrudez(xyz(2), tsuppb, Rbo, Rbi, minsuppheight)/supportlayers)], true );
            
            suppfens = translate(suppfens, [0,0,extrudedir*tbi/2]);

            [fens, gcells, suppgcells] = merge_meshes(shaftfens, shaftgcells, suppfens, suppgcells, mergetol);

            gcells = cat(gcells, suppgcells);

        end
        
        % subdivide mesh into 20 nodes instead of 8 nodes
%         [fens,gcells] = H8_to_H20(fens,gcells);
        
%         fens = shaftfens;
%         gcells = shaftgcells;

        % get the locations of the nodes in the mesh
        xyz = get (fens,'xyz');
        
        [x,y] = pol2cart(-xyz(:,1),xyz(:,2));
    
        xyz = [x, y, xyz(:,3)];
        
%         % shift to match rotor disc positon which is located at the origin
%         xyz(:,3) = xyz(:,3) - (tbi/2);
        
        % replace the coordinates of the nodes
        fens = set(fens, 'xyz', xyz);
        
        % Merge together the nodes at the fin edges joint
        [fens, gcells] = merge_nodes(fens, gcells, mergetol);
        
        % merge with the fin
        [fens, gcells, fingcells] = merge_meshes(fens, gcells, finfens, fingcells, mergetol);
        
        gcells = cat(gcells, fingcells);
        
    else
        
        fens = finfens;
        gcells = fingcells;
        
    end
    
%     drawmesh({fens, gcells});

end


function z = tsuppbextrudez(y, maxh, Rbo, Rbi, tol)
% calculates a appropriate extrusion height to create supports for the
% rotor which reduce in height from the specified maximum at Rbi to zero at
% Rbi

    z = max(tol, -(maxh ./ (Rbo - Rbi)).*y + (maxh.*Rbo ./ (Rbo - Rbi)));

%     z = maxh;
    
end