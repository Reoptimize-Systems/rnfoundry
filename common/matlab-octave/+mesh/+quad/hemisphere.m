function [vertices, faces, info, hfig, hax] = hemisphere (R, n_patch, varargin)
% creates a quad shell mesh of a hemisphere that can optionally be sectioned
%
% Syntax
%
% [vertices, faces, info, hfig, hax] = mesh.quad.hemisphere (R, n_patch)
% [vertices, faces, info, hfig, hax] = mesh.quad.hemisphere (..., 'Parameter', Value)
%
% Description
%
% mesh.quad.hemisphere creates a quadrilateral mesh of a hemisphere with
% the centre (if it where a full sphere) at the origin and the flat edge
% lying on the Z-Y plane. Optionally the hemisphere can also be truncated
% at a desired height in the Z direction.
%
% Input
%
%  R - radius of the hemisphere
%
%  n_patch - number of quads in the patch making up the end of each
%   hemishpere, and determining the overall density of the mesh.
%
% Addtional arguments may be supplied as parameter-value pairs. The
% available options are:
%
%  'CutZMax' - Cut off any faces which extend higher in Z direction than 
%    the value specified here. A new flat edge will be created at this Z
%    position on a plane parallel to the X-Y plane.
%
%  'Draw' - optional true/false flag indicating whether to plot the mesh in
%    a figure. Default is false if not supplied.
%
%  'NJointCircles' - optional scalar integer which is the number of
%    vertical slices making up the mesh of the hemisphere between the end
%    patch and the start of the cylindrical part. The same mesh spacing
%    will then also be used for the cylinder. Default is sqrt(n_patch) if
%    not supplied.
%
% Output
%
%  vertices - (3 x p) matrix of x, y and z coordinated of mesh vertices
%
%  faces - (4 x q) matrix of indices into the vertices matrix. Each colum
%   representing on quadrilateral element of the mesh.
%
%  info - structure containing information about the mesh, the following
%   fields can be present:
%
%   EndInds : indices of the columns of the vertices matrix containing the
%    vertices on the flat face of the sphere on the Z-Y plane. Indices are
%    provided going clockwise from the highest Z point when looking down
%    the X axis in the -ve X direction.
%
%   EndThetas : the angular positions of the vertices on the flat face of
%    the sphere on the Z-Y plane which are referred to in the EndInds
%    field above. The angles are the angular positions of the points on a
%    circle on the Z-Y plane with centre at (0,0,0), with angle of zero
%    being at the point of maximum Z.
%
%   CutTopInds : indices of the columns of the vertices matrix containing 
%    the vertices on the flat face of the sphere when it has been cut
%    through the X-Y plane using the CutZMax option. Indices are provided
%    going clockwise around the new top cut edge from the point of minimum
%    Y to maximum Y.
%
%  hfig - if 'Draw' is true this will contain a handle to the figure
%   containing the mesh plot. If draw is false, it will be the empty matrix
%   [].
%
%  hax -  if 'Draw' is true this will contain a handle to the axes
%   containing the mesh plot. If draw is false, it will be the empty matrix
%   [].
%
%
%
% See Also: mesh.quad.spherepatchgrid
%
    options.CutZMax = [];
    options.CutZMin = [];
    options.CutYMax = [];
    options.CutYMin = [];
    options.NJointCircles = ceil (sqrt(n_patch));
    options.Draw = false;
    
    options = parse_pv_pairs (options, varargin);
    
    hfig = []; 
    hax = [];

    [vertices, faces, sphere_patch_info] = mesh.quad.spherepatchgrid (R, n_patch, ...
                                                    'ProjectionType', 'equidistance', ...
                                                    'CutZMax', options.CutZMax, ...
                                                    'CutZMin', options.CutZMin, ...
                                                    'CutYMax', options.CutYMax, ...
                                                    'CutYMin', options.CutYMin );

    if isempty (vertices)
        error ('No faces were left on the end patch acter cutting, unable to continue');
    end

    if options.Draw
        
        hfig = figure;
        
        hax = axes ('Parent',hfig);

        patch('vertices', vertices', 'faces', faces', 'FaceAlpha', 0, 'parent',hax);

        axis(hax,'equal')
        
    end

    npoints = numel (sphere_patch_info.AllEdgeInds) - numel (sphere_patch_info.TopEdgeInds) + 2;

    
    circ_conn_inds = [ fliplr(sphere_patch_info.LeftEdgeInds(2:end)'), ...
                       sphere_patch_info.BottomEdgeInds', ...
                       sphere_patch_info.RightEdgeInds(2:end)' ];
                   
    if options.CutZMax >= R
        % if the hemisphere is cut off above the top of the end patch, we
        % add in the patch top edge indices
        circ_conn_inds = [ circ_conn_inds, sphere_patch_info.TopEdgeInds(2:end-1)' ];
    end
    
    % we'll now make up the rest of the hemisphere by joining circular
    % slices made parallel to the Z-Y plane
    
    % first find points at the same angular positions as the end patch edge
    % vertices
    end_thetas = cart2pol ( vertices( 2, circ_conn_inds ), vertices( 3, circ_conn_inds) );

    % find the minimum distance of the patch vertices from the Y-Z plane,
    % we will split the distance between this x position and the plane up
    % into the desired number of slices
    lowest_x = min (vertices(1,:));

    % ensure we don't overlap with the patch
    first_circ_pos_x = 0.9 * lowest_x;

    % next find the radius of the first slice of the sphere at the desired
    % X position.
    circ_slice_R = sqrt (R^2 - first_circ_pos_x^2);

    [ circ_conn_Y, circ_conn_Z ] = pol2cart (end_thetas, circ_slice_R);

    new_vertices = [ repmat(first_circ_pos_x, size (circ_conn_Z)); circ_conn_Y; circ_conn_Z ];

    nverts = size (vertices, 2);

    vertices = [ vertices, new_vertices ];

    new_vert_inds = (nverts+1):(size (vertices, 2));

    new_faces = [ circ_conn_inds(1:end-1);
                  new_vert_inds(1:end-1);
                  new_vert_inds(2:end);
                  circ_conn_inds(2:end) ]; 
              
	faces = [ faces, new_faces ];

    top_inds = [ new_vert_inds(1), new_vert_inds(end) ];

    if options.Draw
        hold on;
        patch('vertices', vertices', 'faces', new_faces', 'FaceAlpha', 0, 'EdgeColor', 'b', 'parent',hax);
        hold off;
    end

    prev_vert_inds = new_vert_inds;

    for ind = options.NJointCircles-1:-1:0

        circ_pos_x = (ind/options.NJointCircles) * first_circ_pos_x;

        circ_slice_R = sqrt (R^2 - circ_pos_x^2);

        [ circ_conn_Y, circ_conn_Z ] = pol2cart (end_thetas, circ_slice_R);

        new_vertices = [ repmat(circ_pos_x, size (circ_conn_Z)); circ_conn_Y; circ_conn_Z ];

        nverts = size (vertices, 2);

        vertices = [ vertices, new_vertices ];

        new_vert_inds = (nverts+1):(size (vertices, 2));

        new_faces = [ prev_vert_inds(1:end-1);
                      new_vert_inds(1:end-1);
                      new_vert_inds(2:end);
                      prev_vert_inds(2:end) ];
                  
        faces = [ faces, new_faces ];

        if options.Draw
            hold on;
            patch('vertices', vertices', 'faces', new_faces', 'FaceAlpha', 0, 'EdgeColor', 'b', 'parent',hax);
            hold off;
        end

        top_inds = [ top_inds; new_vert_inds(1), new_vert_inds(end) ];

        prev_vert_inds = new_vert_inds;

    end
    
    info.EndInds = new_vert_inds;

    top_inds = flipud (top_inds);

    if sphere_patch_info.FullPatch
        all_top_inds = [ top_inds(:,1)', sphere_patch_info.TopEdgeInds', flipud(top_inds(:,2))' ];
    else
        all_top_inds = [ top_inds(:,1)', sphere_patch_info.TopEdgeInds', flipud(top_inds(:,2))' ];
    end

    % now draw the top edge and faces
    
    % get Y position of cut at the X axis to give the radius of the slice
    % through the sphere at this point
    circ_slice_R = sqrt (R^2 - (options.CutZMax)^2);

    top_thetas = cart2pol ( vertices( 1, all_top_inds ), vertices( 2, all_top_inds) );

    [ circ_conn_X, circ_conn_Y ] = pol2cart (top_thetas, circ_slice_R);

    if options.Draw
        hold on
        scatter3 (circ_conn_X, circ_conn_Y, repmat(options.CutZMax, size (circ_conn_X)), 'ob');
        hold off
    end

    new_vertices = [ circ_conn_X; circ_conn_Y; repmat(options.CutZMax, size (circ_conn_X)) ];

    nverts = size (vertices, 2);

    vertices = [ vertices, new_vertices ];

    info.CutTopInds = (nverts+1):(size (vertices, 2));

    info.EndInds = [ info.CutTopInds(1), info.EndInds, info.CutTopInds(end) ];
    
    new_faces = [ all_top_inds(2:end);
                  info.CutTopInds(2:end);
                  info.CutTopInds(1:end-1);
                  all_top_inds(1:end-1);
                ];

	faces = [ faces, new_faces ];
    
    info.EndThetas = cart2pol ( vertices( 2, info.EndInds ), vertices( 3, info.EndInds) );
    
    if options.Draw
        hold on;
        patch('vertices', vertices', 'faces', new_faces', 'FaceAlpha', 0, 'EdgeColor', 'g', 'parent',hax);
        hold off;
    end


end