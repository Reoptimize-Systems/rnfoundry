function [fens,gcells] = axialrotorfinlinkmesh3d(Rs, Rbi, tbi, outersep, modulexs, nmodules, nlinks, layers, varargin)
% creates a structured FAESOR mesh of the shaft portions between fins in an
% axial flux rotor
%
% Syntax
% 
% [fens,gcells] = axialrotorfinlinkmesh3d(Rs, Rbi, tbi, outersep, modulexs, nlinks, layers)
%
% Input
%
% Rs - inner radius of rotor shaft 
% 
% Rbi - outer radius of rotor shaft, and inner radius of the rotor fin back
%   iron
% 
% tlink - axial length of the link
% 
% circumpoints - number of circumferential points the mesh will be made up
%   of (before refinement)
%
% cellsperm - number cells per metre of axail length in the link
%
% Output
%
% fens - FAESOR fenodeset object containing the nodes of the shaft links mesh
%
% gcells - FAESOR gcellset_H20 oject containing the H20 type cells making
%   up the rotor shaft links mesh
%
%

    Inputs.ShaftLabels = (10:10+nlinks);
    Inputs.Subset = nmodules;
    if nmodules == 1
        Inputs.ModuleGap = 0;
    else
        Inputs.ModuleGap = 2*pi/(100*nmodules);
    end
    
    Inputs = parseoptions(Inputs, varargin);
    
    if isscalar(Inputs.ShaftLabels)
        Inputs.ShaftLabels = repmat(Inputs.ShaftLabels, 1, nlinks);
    elseif numel(Inputs.ShaftLabels) ~= nlinks
        error('Number of shaft labels does not match number of links and is not a scalar.')
    end
    
    tshaftsec = (outersep - tbi*(nlinks - 1)) / (nlinks);
    
    % make a base link
    [basefens,basegcells] = axialrotorshaftmesh3d(Rs, Rbi, tshaftsec, ...
            nmodules, layers, modulexs, 'Subset', Inputs.Subset, 'ModuleGap', Inputs.ModuleGap);
       
    basegcells = setlabels(basegcells, Inputs.ShaftLabels(1));
    
    fens = basefens;
    gcells = basegcells;
    
    for i = 2:nlinks
        
        % copy over the base gcells
        nextfens = basefens;
        nextgcells = basegcells;
        
        % shift the cells up
        xyz = get(nextfens, 'xyz');
        xyz(:,3) = xyz(:,3) + (i-1)*(tshaftsec + tbi);
        nextfens = set(nextfens, 'xyz', xyz);
        
        nextgcells = setlabels(nextgcells, Inputs.ShaftLabels(i));
        
        % merge the meshes
        % merge the two meshes into one
        [fens,gcells,nextgcells] = merge_meshes(fens, gcells, nextfens, nextgcells, 0);

        % concatenate the cells into one for convenience
        gcells = cat(gcells, nextgcells);
        
    end
    
    % shft every think to the middle
    xyz = get(fens, 'xyz');
    xyz(:,3) = xyz(:,3) - outersep/2 + tshaftsec/2;
    fens = set(fens, 'xyz', xyz);
    
end