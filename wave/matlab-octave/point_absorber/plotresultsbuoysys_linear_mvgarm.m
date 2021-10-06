function plotresultsbuoysys_linear_mvgarm(T, Y, results, design, skip)
% Syntax:
%
%   plotresultsbuoysys_linear_mvgarm(T, Y, results, design, skip)
%
% Input:
%
%   T - vector of time steps for which the values of the the snapper
%       outputs in Y have been evaluated
%
%   Y - Matrix containing the the buoy position in heave, the buoy velocity
%       in heave, the buoy position in surge, the buoy velocity in surge,
%       the armature displacment, velocity and coil current, and at each
%       time step in the vector T. So that:
%
%       Y = [ xBh, vBh, xBs, vBs, I1, I2, I3, ...]
%
%   results - A structure containing some further outputs from the
%             simulation. It should contain a the following fields: 
%
%             EMF: a vector of the EMF values at each time step in the
%             vector T. 
%
%             Ffea: the electromagnetic forces at each time step in T
%
%             Fs: the spring forces at each time step in T
%
%             ydot: an (n x p) matrix containing the results of the 
%             right hand side of the differential equations solved by
%             the matlab ode solvers.
%
%             xT: Translator position
%
%             vT: Translator velocity
%
%             Ff: (optional) The frictional forces between the armature and
%             translator at each time step in T.
%
%   skip - A value determining what proportion of time points are to be
%          plotted. If skip is 2, every other step is plotted, if 3, every
%          third point etc. 


    if nargin < 5
        skip = 1;
    end
    
    % Y(:,1) is the buoy displacement in heave
    % Y(:,2) is the buoy velocity in heave
    % Y(:,3) is the buoy displacement in surge
    % Y(:,4) is the buoy velocity in surge    
    % Y(:,5) is the field displacement
    % Y(:,6) is the field velocity
    % Y(:,7) is the first coil current
    xBh = Y(:,1);
    vBh = Y(:,2);
    xBs = Y(:,3);
    vBs = Y(:,4);
    xA = Y(:,5);
    vA = Y(:,6);
    Icoils = Y(:,7:6+design.phases);

    % create the figure, ensuring it is visible
%     h = figure('visible','on', 'Units','normalized','outerposition',[0 0.02 1 0.98]);
    
    h = figure('visible','on', 'Units','normalized');
    
    if isoctave
        set(h, 'Position', [0, 0.05, 0.95, 0.95]);
    else
        maximize(h);
    end

    % we will create a 2 x 1 subplot
    hax = subplot(2,1,1);
    
    [hax, legendstrings] = plotmachineresults_linear(T, 'xA', xA, ...
                                                        'vA', vA, ...
                                                        'Icoils', Icoils, ...
                                                        'EMF', results.EMF, ...
                                                        'vT', results.vT, ...
                                                        'xT', results.xT, ...
                                                        'AxesPosition', [0.14/2, 0.55, 0.86, 0.42], ...
                                                        'hax', hax);
    

    % Create a plot legend
    lh = legend(legendstrings{:}, 'Location', 'NorthWest');
    
    set(lh, 'Box', 'off');
    set(lh, 'FontSize', 10);

    pos = get(hax, 'Position');

    subplot(2,1,2);

    legendstrings = plotforces_linear(T, results, skip);
    
    lh = legend(legendstrings, 'Location', 'NorthWest');
    set(lh, 'Box', 'off');%legend(lh, 'boxoff');
    set(lh, 'FontSize', 10);
    
    hax = gca;
    set(hax, 'Position', [pos(1), 0.05, pos(3), 0.42])
    
    set(lh, 'Box', 'off');%legend(lh, 'boxoff');
    set(lh, 'FontSize', 10);

    % Now plot some hydrodynamic results
    plotbuoyresults_linear(T, xBh, vBh, xBs, vBs, results, skip);
    
end