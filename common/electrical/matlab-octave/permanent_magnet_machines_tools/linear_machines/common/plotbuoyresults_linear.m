function hfigs = plotbuoyresults_linear(T, xBh, vBh, xBs, vBs, results, skip)
% plots some results from the hydrodynamic model which are of interest
%
% Syntax:
%
%       plotbuoyresults_linear(T, xBh, vBh, xBs, vBs, results, skip)
%
% Input:
%
%
% Output:
%
%
 
    % Create a new figure to show positions
    hfigs = figure;
    set(hfigs(1), 'Units', 'Normalized', 'Position', [0, 0.05, 0.85, 0.85]);
    
    % plot the wave height
    plot (T, results.wave_height)
    
    % Add new plots to the same graphs
    hold on
    
    % plot the buoy position in heave
    plot(T, xBh, '-r')
    
    % plot the buoy position in surge
    plot(T, xBs, ':k')
    
    % plot the translator position
    if isfield(results, 'xE')
        plot(T, results.xE, ':gx', 'MarkerSize', 2)
    elseif isfield(results, 'xT')
        plot(T, results.xT, ':gx', 'MarkerSize', 2)
    end
    
    hold off
    
    % Set the graph titles etc.
    title('Positions')
    xlabel('Time (s)');
    ylabel('Position (m)');
    legend('Wave Height', 'Buoy Position -h', 'Buoy Position -s', ...
        'Translator Position -h')
    
    % Create a new figure to show velocities
    hfigs(2) = figure;
    set(hfigs(2), 'Units', 'Normalized', 'Position', [0, 0.05, 0.85, 0.85]);
    
    % Add new plots to the same graphs
    hold on
    
    % plot the buoy position in heave
    plot(T, vBh, '-r')
    
    % plot the buoy position in surge
    plot(T, vBs, ':k')
    
    % plot the translator position
    if isfield(results, 'vE')
        plot(T, results.vE, ':gx', 'MarkerSize', 2)
    elseif isfield(results, 'vT')
        plot(T, results.vT, ':gx', 'MarkerSize', 2)
    end
    
    hold off
    
    % Set the graph titles etc.
    title('Velocities')
    xlabel('Time (s)');
    ylabel('Velocity (ms^{-1})');
    legend('Buoy Velocity -h', 'Buoy Velocity -s', 'Translator Velocity -h')

    % Plot new figure showing the forces on the bodies
    hfigs(3) = figure;
    set(hfigs(3), 'Units', 'Normalized', 'Position', [0, 0.05, 0.85, 0.85]);

    legstr = {};
    
    plot(T(1:skip:length(T)), results.buoyancy_force(1:skip:length(T)))
    
    legstr = [legstr, {'Buoyancy'}];
    
    hold on
    
    plot(T(1:skip:length(T)), results.FBDh(1:skip:length(T)),'m')
    
    legstr = [legstr, {'Buoy drag force -h'}];
    
    plot(T(1:skip:length(T)), results.FBDs(1:skip:length(T)),':m')
    
    legstr = [legstr, {'Buoy drag force -s'}];
    
    plot(T(1:skip:length(T)), results.excitation_force_heave(1:skip:length(T)), 'r')
    
    legstr = [legstr, {'Excitation -h'}];
    
    plot(T(1:skip:length(T)), results.excitation_force_surge(1:skip:length(T)), ':r')
    
    legstr = [legstr, {'Excitation -s'}];
    
    plot(T(1:skip:length(T)), results.radiation_force_heave(1:skip:length(T)), 'b')
    
    legstr = [legstr, {'Radiation -h'}];
    
    plot(T(1:skip:length(T)), results.radiation_force_surge(1:skip:length(T)), ':b')
    
    legstr = [legstr, {'Radiation -s'}];
    
    netforce_heave = results.excitation_force_heave(1:skip:length(T)) + results.buoyancy_force(1:skip:length(T)) ...
                        + results.radiation_force_heave(1:skip:length(T)) + results.FBDh(1:skip:length(T));
                    
    netforce_surge = results.excitation_force_surge(1:skip:length(T)) + results.radiation_force_surge(1:skip:length(T)) + results.FBDs(1:skip:length(T));              
    
    if isfield(results, 'Ffea_surge') && isfield(results, 'Ffea_heave')
        plot(T(1:skip:length(T)), results.Ffea_heave(1:skip:length(T)), 'g')
        plot(T(1:skip:length(T)), results.Ffea_surge(1:skip:length(T)), ':g')
        legstr = [legstr, {'Ffea -h', 'Ffea -s'}];
        netforce_heave = netforce_heave + results.Ffea_heave(1:skip:length(T));
        netforce_surge = netforce_surge + results.Ffea_surge(1:skip:length(T));
    end
    
    if isfield(results, 'FaddB')
        plot(T(1:skip:length(T)), results.FaddB(1:skip:length(T),1), 'y')
        plot(T(1:skip:length(T)), results.FaddB(1:skip:length(T),2), ':y')
        legstr = [legstr, {'FaddB -h', 'FaddB -s'}];
        netforce_heave = netforce_heave + results.FaddB(1:skip:length(T),1);
        netforce_surge = netforce_surge + results.FaddB(1:skip:length(T),2);
    end

    % Net force on buoy in heave
    plot(T(1:skip:length(T)), netforce_heave, 'k');
    
    legstr = [legstr, {'Net Forces -h'}];
    
    % Net force on buoy in surge
    plot(T(1:skip:length(T)), netforce_surge, ':k')
    
    legstr = [legstr, {'Net Forces -s'}];

    hold off
    
    title('System Forces')
    xlabel('Time (s)');
    ylabel('Forces (N)');
    
    legend(legstr{:});
           
    
end