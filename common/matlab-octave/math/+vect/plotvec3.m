function [hquiv, hax] = plotvec3 (vec, startpoints, varargin)
    % plot one or more 3D vectors
    %
    % Syntax
    %
    % plotvec3 (vec, startpoints, 'Parameter', value, ...)
    %
    % Input
    %
    %  vec - (3 x n) matrix of one or more 3D vectors to be plotted
    %
    %  startpoints - optional (3 x n) matrix of start points for the
    %    vectors in vec, arrows of the sizes in vec will be drawn from
    %    startpoints. If not supplied or empty, all vectors in vec will be
    %    drawn from the origin.
    %
    %  Additional arguments may be supplied as parameter-value pairs. The
    %  available options are:
    %
    %  'PlotAxes' - handle to axes in which to do the plot, if not
    %    supplied, a new set of axes is created;
    %
    %  'Scale' - automatically scales the vectors to prevent them from
    %    overlapping, and then multiplies them by scale. scale = 2 doubles
    %    their relative length, and scale = 0.5 halves them. Use scale = 0
    %    to plot the vectors without the automatic scaling. Default if not
    %    supplied is 1.
    %
    %  'Properties' - cell array of properties to be passed to quiver3
    %    which is used to make the plot. This should be a set of
    %    parameter-value pairs
    %
    %

    options.PlotAxes = [];
    options.Scale = 1;
    options.Properties = {};
    
    options = parse_pv_pairs (options, varargin);
    
    if nargin < 2 || isempty (startpoints)
        startpoints = repmat (zeros, size(vec));
    end
    
    if isempty (options.PlotAxes)
        hax = axes;
    else
        hax = options.PlotAxes;
    end
    
    hold on
    
    hquiv = quiver3 ( hax, startpoints(1,:), ...
                      startpoints(2,:), ...
                      startpoints(3,:), ...
                      vec(1,:), ...
                      vec(2,:), ...
                      vec(3,:), ...
                      options.Scale, ...
                      options.Properties{:} );

	hold off

end