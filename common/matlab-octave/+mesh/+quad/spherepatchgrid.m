function [vertices, faces, info] = spherepatchgrid (R, n, varargin)
% creates a quad mesh of a patch on a sphere by projecting a grid onto it
%
% Syntax
%
% [vertices, faces, edge_inds] = cubedspherepatch (R, n)
% [vertices, faces, edge_inds] = cubedspherepatch (..., 'Parameter', Value)
%
% Description
%
% cubedspherepatch creates a quad mesh of the end of a sphere which can
% optionally be clipped above a certain Z or Y position. It works by
% projecting a grid onto the surface of a shere, and is useful for creating
% quadrilateral meshes of spherically curved domains, without degenerate
% quads (i.e triangles).
%
% Input
%
%  R - radius of the sphere end to be meshed
%
%  n - number of faces to produce on the end of the sphere
%
% Addtional arguments may be supplied as parameter-value pairs. The
% available options are:
%
%  'ProjectionType' - character vector indicating the type of projection to
%    use
%
%  'CutZMax' - Cut off any faces which extend higher in Z direction than 
%    the value specified here
%
%  'CutZMin' - Cut off any faces which extend lower in Z direction than 
%    the value specified here
%
%  'CutYMax' - Cut off any faces which extend in positive Y direction more  
%    than the value specified here
%
%  'CutYMin' - Cut off any faces which have vertices with Y positions less  
%    than the value specified here
%
% Output
%
%  vertices - (3 x v) matrix of vertices
%
%  faces - (4 x u) matrix of indices into the vertices matrix where each
%    column represents one face
%
%  info - structure contating further information on the mesh. It will 
%    contain the following fields:
%
%    AllEdgeInds : vector of indices into the vertices matrix for the 
%      points on the outer edge of the generated mesh, these are the points
%      on the outer boundary starting from the the point with the minimum Y
%      position and minimum Z position and moving anti-clockwise around the
%      outer boundary.
%
%    BottomEdgeInds - vector of indices into the vertices matrix for the 
%      points on the lower outer edge of the generated mesh, i.e. the edge
%      points at the minimu Z position. Arranged from minimum Y position to
%      maximum Y position.
%
%    RightEdgeInds - vector of indices into the vertices matrix for the 
%      points on the right outer edge of the generated mesh, i.e. the edge
%      points at the maximum Y position. Arranged from minimum Z position
%      to maximum Z position.
%
%    TopEdgeInds - vector of indices into the vertices matrix for the 
%      points on the upper outer edge of the generated mesh, i.e. the edge
%      points at the maximum Z position. Arranged from maximum Y position
%      to minimum Y position.
%
%    LeftEdgeInds - vector of indices into the vertices matrix for the 
%      points on the lower outer edge of the generated mesh, i.e. the edge
%      points at the minimum Y position. Arranged from maximum Z position
%      to minimum Z position.
%


    options.CutZMax = [];
    options.CutZMin = [];
    options.CutYMax = [];
    options.CutYMin = [];
    
    options.ProjectionType = 'equidistance';
    
    options = parse_pv_pairs (options, varargin);
    
    n = round (n);
    
    if n < 1
        error ('cubedsphere: n must be stricly positive number');
    end
    
    n = n+1;
    
    % Discretize basic face of a cube
    switch lower (options.ProjectionType)
        case 'equiangular'
            % equidistance projection
            x = R;
            theta = linspace (-pi/4, pi/4, n);
            y = tan (theta);
            z = y;
            [X, Y, Z] = ndgrid (x,y,z);
        case 'equidistance'
            % equidistance projection
            x = R;
            y = linspace (-R, R, n);
            z = y;
            [X, Y, Z] = ndgrid (x,y,z);
        otherwise
            error('cubedsphere: ProjectionType must be ''equiangular'' or ''equidistance''');
    end
    
    % Patch topology
    [i, j] = meshgrid (1:n-1,1:n-1); % 090815: Change to meshgrid from ndgrid
    i = i(:); j=j(:);
    faces = [ sub2ind([n n], i, j) ...
              sub2ind([n n], i+1, j) ...
              sub2ind([n n], i+1, j+1) ...
              sub2ind([n n], i, j+1) ].';
         

    X = X(:);
    Y = Y(:);
    Z = Z(:);
    
    S = R*sqrt(1./(X.^2+Y.^2+Z.^2));
    
%     X_tmp = X.*S;
    Y_tmp = Y.*S;
    Z_tmp = Z.*S;
    
%     remove_face = [];
%     uncut_inds = [];

    info.FullPatch = true;
    
    if ~isempty (options.CutZMax)
        
        % gow down the rows from the top looking for rows with any 
        % Z >= CutZMax
        for rowind = n:-1:1
            
            this_row_vert_inds = (rowind-1)*n+1:rowind*n;
            
            if any (Z_tmp(this_row_vert_inds) >= options.CutZMax)
                
                info.FullPatch = false;

                % remove the whole row of points and faces
                cut_inds = this_row_vert_inds;

                [X, Y, Z, faces ] = remove_faces_and_inds_and_renumber ( X, Y, Z, faces, cut_inds );
                
%                 remove_face = [];
%                 for face_ind = 1:size(faces, 2)
%                     for vert_ind = 1:4
%                         if any (faces(vert_ind,face_ind) == cut_inds)
%                             remove_face = [remove_face, face_ind];
%                         end
%                     end
%                 end
% 
%                 faces(:,remove_face) = [];
% 
%                 X(cut_inds) = [];
%                 Y(cut_inds) = [];
%                 Z(cut_inds) = [];
            
            else
                % none of the points on this row are above the threshold,
                % so none below will be either, so stop searching
                break;
            end
        
        end
        
    end
    
    if ~isempty (options.CutZMin)
       
        uncut_inds = [ uncut_inds, find(Z_tmp >= options.CutZMin)];
        
        cut_inds = find (Z_tmp < options.CutZMin);
        
        for face_ind = 1:size(faces, 2)
            for vert_ind = 1:4
                if any (faces(vert_ind,face_ind) == cut_inds)
                    remove_face = [remove_face, face_ind];
                end
            end
        end
        
%         faces(:,remove_face) = [];
%         
%         X = X(uncut_inds);
%         Y_tmp = Y_tmp(uncut_inds);
%         Z_tmp = Z_tmp(uncut_inds);
        
    end
    
    if ~isempty (options.CutYMax)
       
        uncut_inds = [ uncut_inds, find(Y_tmp <= options.CutYMax) ];
        
        cut_inds = find (Y_tmp > options.CutYMax);
        
        for face_ind = 1:size(faces, 2)
            for vert_ind = 1:4
                if any (faces(vert_ind,face_ind) == cut_inds)
                    remove_face = [remove_face, face_ind];
                end
            end
        end
        
%         faces(:,remove_face) = [];
%         
%         X = X(uncut_inds);
%         Y_tmp = Y_tmp(uncut_inds);
%         Z_tmp = Z_tmp(uncut_inds);
        
    end
    
    if ~isempty (options.CutYMin)
       
        uncut_inds = [ uncut_inds, find(Y_tmp >= options.CutYMin) ];
        
        cut_inds = find (Y_tmp < options.CutYMin);
        
        for face_ind = 1:size(faces, 2)
            for vert_ind = 1:4
                if any (faces(vert_ind,face_ind) == cut_inds)
                    remove_face = [remove_face, face_ind];
                end
            end
        end
        
%         faces(:,remove_face) = [];
%         
%         X = X(uncut_inds);
%         Y_tmp = Y_tmp(uncut_inds);
%         Z_tmp = Z_tmp(uncut_inds);
        
    end

%     if ~isempty (remove_face)
%         faces(:,remove_face) = [];
%     end
% 
%     if ~isempty (uncut_inds)
%         X = X(uncut_inds);
%         Y = Y(uncut_inds);
%         Z = Z(uncut_inds);
%     end
    
    info.AllEdgeInds = [];
    
    info.BottomEdgeInds = find (Z <= min(Z));

    tmpY2 = Y(info.BottomEdgeInds);

    [~, I] = sort (tmpY2, 'ascend');
    
    info.AllEdgeInds = [ info.AllEdgeInds; info.BottomEdgeInds(I)];
    
    info.RightEdgeInds = find (Y >= max(Y));
    
    tmpZ2 = Z(info.RightEdgeInds);
    
    [~, I] = sort (tmpZ2, 'ascend');
    
    info.AllEdgeInds = [ info.AllEdgeInds; info.RightEdgeInds(I(2:end))];
    
    info.TopEdgeInds = find (Z >= max(Z));
    
    tmpY2 = Y(info.TopEdgeInds);
    
    [~, I] = sort (tmpY2, 'descend');
    
    info.AllEdgeInds = [ info.AllEdgeInds; info.TopEdgeInds(I(2:end))];
    
    info.LeftEdgeInds = find (Y <= min(Y));
    
    tmpZ2 = Z(info.LeftEdgeInds);
    
    [~, I] = sort (tmpZ2, 'descend');
    
    info.AllEdgeInds = [ info.AllEdgeInds; info.LeftEdgeInds(I(2:end-1))];


    % Project on sphere S2
    S = R*sqrt(1./(X.^2+Y.^2+Z.^2));
    X = X.*S;
    Y = Y.*S;
    Z = Z.*S;

    vertices = [ X(:), Y(:), Z(:) ].';

end


function [X, Y, Z, faces ] = remove_faces_and_inds_and_renumber ( X, Y, Z, faces, cut_inds )

    remove_face = [];
    for face_ind = 1:size(faces, 2)
        for vert_ind = 1:4
            if any (faces(vert_ind,face_ind) == cut_inds)
                remove_face = [remove_face, face_ind];
            end
        end
    end

    faces(:,remove_face) = [];
    
    cut_inds = sort (cut_inds, 'descend');
    
    for cut_inds_ind = 1:numel (cut_inds)
        
        faces (faces > cut_inds (cut_inds_ind)) = faces (faces > cut_inds (cut_inds_ind)) - 1;

        X(cut_inds(cut_inds_ind)) = [];
        Y(cut_inds(cut_inds_ind)) = [];
        Z(cut_inds(cut_inds_ind)) = [];
    
    end

end