function [fens, gcells] = axialinnerrotorfinsmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, outersep, nmodules, nfins, maglabels, layers, modulexs, varargin)
% creates a structured FAESOR mesh of the inner rotor fins of a of an
% axial-flux type machine
%
% Syntax
% 
% [fens,gcells] = axialinnerrotorfinsmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, layers, circumpoints)
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
% separation - this is the space between the inner sufaces of the outer
%   fins
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

    Options.BackIronRadialPointsPerM = 15;
    Options.MagnetRadialPointsPerM = 15;
    Options.Subset = nmodules;
    if nmodules == 1
        Options.ModuleGap = 0;
    else
        Options.ModuleGap = 2*pi/(100*nmodules);
    end
    
    Options = parseoptions(Options, varargin);
                              
    if nargin < 9 || isempty(maglabels)
        maglabels = zeros(1, nfins);
    end
    
    % get the fin positions
    pos = linspace(0, outersep+tbi, nfins+2);
    
    for i = 2:numel(pos)-1
        
        [newfens,newgcells] = axialrotorfinmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, nmodules, layers, modulexs, maglabels(i-1), ...
                                  'BackIronRadialPointsPerM', Options.BackIronRadialPointsPerM, ...
                                  'MagnetRadialPointsPerM', Options.MagnetRadialPointsPerM, ...
                                  'Subset', Options.Subset, ...
                                  'ModuleGap', Options.ModuleGap);

        % move the mesh in the -ve x direction 
        xyz = get(newfens, 'xyz');
        xyz(:,3) = xyz(:,3) + pos(i);
        newfens = set(newfens, 'xyz', xyz);
    
        if i > 2
            
            % merge the two meshes into one
            [fens,gcells,newgcells] = merge_meshes(fens, gcells, newfens, newgcells, 0);
            
            % concatenate the cells
            gcells = cat(gcells, newgcells);
            
        else
            fens = newfens;
            gcells = newgcells;
        end
    
    end

    % move the mesh to the middle
    xyz = get(fens, 'xyz');
    xyz(:,3) = xyz(:,3) - ((outersep+tbi)/2);
    fens = set(fens, 'xyz', xyz);
    
end