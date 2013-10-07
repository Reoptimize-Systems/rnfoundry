function plotresultsproscribedmot_linear_mvgarm(T, Y, results, design, skip)
%   T - vector of time steps for which the values of the the snapper
%       outputs in Y have been evaluated
%
%  Y - This is the output from the ode solver for a proscribed motion
%      simulation. In this case the vector contains the armature
%      displacment, velocity and coil current respectively at each time
%      step in the vector T. So that:
%
%       Y = [ Iphases, xA, vA ]
%
%   results - A structure containing some further outputs from the
%             simulation. It should contain the following fields: 
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
%             The following fields can also optionally be supplied:
%
%             FDrive: The drive force applied to the translator at each
%             time step in T.
%
%             xT: Translator position
%
%             vT: Translator velocity
%
%             Ff: (optional) The frictional forces between the armature and
%             translator at each time step in T.
%

%
%   skip - A value determining what proportion of time points are to be
%          plotted. If skip is 2, every other step is plotted, if 3, every
%          third point etc. 

    % Y is the coil current for each phase, where each phase has it's own
    % column
    
    Icoils = Y(:,1:design.Phases);
    
    xF = Y(:,design.Phases + 1);
    vF = Y(:,design.Phases + 2);
    
%     
%     % Get the Y-axis limits for the various plots
%     maxY = [1.05*max([abs([Icoils, results.xT, results.vT, xA, vA]); repmat(1e-10,1,size(Icoils,2)+4)], [], 1)];

    % create the figure, ensuring it is visible
    h = figure('visible','on', 'Units','normalized','outerposition',[0 0.02 1 0.98]);
    
    % we will create a 2 x 1 subplot
    hax = subplot(2,1,1);
    
    
    
%     % First plot the positions on the top subplot
%     plot(T(:), results.xT, '--r');
%     hax = gca;
%     set(gca, 'Position', [0.14/2, 0.55, 0.86, 0.42])
%     currentylim = [-1.1*max(results.xT),  1.1*max(results.xT)];
%     ylim(currentylim);
%     ylabel('Displacement, m');%,'FontSize',20);
%     xlabel('Time (s)');%,'FontSize',20);
%     legstr = {'xT : Drive/Translator Position, m'};
% 
%     if ycols > design.Phases
%         
%         hold on
%         
%         plot(T(:), xA, '--b');
% 
%         hold off
%         
%         currentylim = [min(currentylim(1), -1.1*max(abs(xA))), max(currentylim(2), 1.1*max(abs(xA)))];
%         ylim(currentylim);
%         legstr = [legstr, {'xA : Armature Position, m'}];
%         
%     end
%     
%     if ycols > design.Phases
%         
%         currentylim = [-1.1*max(abs(vA)), 1.1*max(abs(vA))];
%         
%     end
% 
%     % Now plot the velocities on their own axes
%     addaxis(T(:), results.vT, [min(currentylim(1), -1.1*max(results.vT)), max(currentylim(2), 1.1*max(results.vT))], '-r');
%     addaxislabel(2,'Velocity, ms^{-1}');
% 
%     if ycols > design.Phases
%        
%         addaxisplot(T(:), vA, 2, '-b');
%         
% %         currentylim = getaddaxisprops(h, 'ylim')
% %         currentylim = [min(currentylim(1), -1.1*max(abs(vA))), max(currentylim(2), 1.1*max(abs(vA)))];
% %         setaddaxisprops(h, 'YLim', currentylim);
% 
%         legstr = [legstr, {'vA : Armature Velovity, m'}];
%         
%     end
%     
%     % Now plot the current on its own axes
%     addaxis(T(:), Icoils(:,1), [-1.1*maxY(1), 1.1*maxY(1)], '-k');
%     addaxislabel(3,'Current, A')
% 
%     % Now plot the EMF on its own axes
%     addaxis(T(1:skip:length(T)), results.EMF(:,1),  [-max(abs([results.EMF(:,1); 1e-10])), max(abs([results.EMF(:,1); 1e-10]))], '-','color', rgb('DarkGreen'));
%     addaxislabel(4,'EMF, V')
% 
% %     lh = legend('xT : Drive/Translator Position, m', 'vT : Drive/Translator Velocity, ms^{-1}', 'I : Coil Current, A', 'EMF, V', 'Location','NorthWest');
%     lh = legend(legstr{:}, 'I : Coil Current, A', 'EMF, V', 'Location','NorthWest');
    
    [hax, legendstrings] = plotmachineresults_linear(T, 'Icoils', Icoils, ...
                                                     'EMF', results.EMF, ...
                                                     'vT', results.vT, ...
                                                     'xT', results.xT, ...
                                                     'xF', xF, ...
                                                     'vF', vF, ...
                                                     'AxesPosition', [0.14/2, 0.55, 0.86, 0.42], ...
                                                     'hax', hax);
                                                    
%     lh = legend(legendstrings, 'Location','NorthWest','orientation','horizontal','-DynamicLegend');
    
    lh = legend(legendstrings, 'Location', 'NorthWest');
    
    set(lh, 'Box', 'off');
    set(lh, 'FontSize', 8);

    pos = get(hax, 'Position');

    % create a subplot for the forces in the results structure
    subplot(2,1,2);

%     % Plot the forces and acceleration on the same graph, starting
%     % with the forces
%     legendstrings = {};
%     currentylim = [0,0];
%     
%     FTtotal = zeros(size(results.Fpto));
%     
%     if isfield(results, 'Fpto')
%         plot(T(1:skip:length(T)), results.Fpto);
%         currentylim = [-1.1*max(abs(results.Fpto)), 1.1*max(abs(results.Fpto))];
%         ylim(currentylim);
%         legendstrings = [legendstrings; {'PTO Force on Translator'}];
%         
%         FTtotal = FTtotal + results.Fpto;
%     end
%     
%     % if present plot the spring forces
%     if isfield(results, 'Fs')
%         hold on
%         plot(T(1:skip:length(T)), results.Fs);
%         hold off
%         currentylim = [min(currentylim(1), -1.1*max(abs(results.Fs))), max(currentylim(2), 1.1*max(abs(results.Fs)))];
%         ylim(currentylim);
%         legendstrings = [legendstrings; {'Spring Force'}];
%     end
%     
%     if isfield(results, 'FfA')
%         hold on
%         plot(T(1:skip:length(T)), results.FfA, 'r');
%         hold off
%         currentylim = [min(currentylim(1), -1.1*max(abs(results.FfA))), max(currentylim(2), 1.1*max(abs(results.FfA)))];
%         ylim(currentylim);
%         legendstrings = [legendstrings; {'Armature Friction'}];
%     end
%     
%     if isfield(results, 'FfT')
%         hold on
%         plot(T(1:skip:length(T)), results.FfT, 'r');
%         hold off
%         currentylim = [min(currentylim(1), -1.1*max(abs(results.FfT))), max(currentylim(2), 1.1*max(abs(results.FfT)))];
%         ylim(currentylim);
%         legendstrings = [legendstrings; {'Armature Friction'}];
%         
%         FTtotal = FTtotal + results.FfT;
%     end
%     
%     if isfield(results, 'Fsnap')
%         hold on
%         plot(T(1:skip:length(T)), results.Fsnap, 'k');
%         hold off
%         currentylim = [min(currentylim(1), -1.1*max(abs(results.Fsnap))), max(currentylim(2), 1.1*max(abs(results.Fsnap)))];
%         ylim(currentylim);
%         legendstrings = [legendstrings; {'Snapping Force (Armature)'}];
%         
%         FTtotal = FTtotal - results.Fsnap;
%         
%     end
%     
%     if isfield(results, 'FLinearDrag')
%         hold on
%         plot(T(1:skip:length(T)), results.FLinearDrag, 'm');
%         hold off
%         currentylim = [min(currentylim(1), -1.1*max(abs(results.FLinearDrag))), max(currentylim(2), 1.1*max(abs(results.FLinearDrag)))];;
%         ylim(currentylim);
%         legendstrings = [legendstrings; {'Linear Drag Force (Armature)'}];
%         
%         FTtotal = FTtotal - results.FLinearDrag;
%     end
%     
%     if isfield(results, 'Faddtrans')
%         
%         hold on
%         plot(T(1:skip:length(T)), results.Faddtrans(:,1), ':k');
%         hold off
%         currentylim = [min(currentylim(1), -1.1*max(abs(results.Faddtrans(:,1)))), max(currentylim(2), 1.1*max(abs(results.Faddtrans(:,1))))];
%         ylim(currentylim);
%         legendstrings = [legendstrings; {'Total Translator Force'}];
%     end
% 
%     xlabel('Time (s)');%,'FontSize',20);
%     ylabel('Force (N)');%, 'FontSize', 20);

    legendstrings = plotforces_linear(T, results, skip);
    
    lh = legend(legendstrings, 'Location', 'NorthWest');
    
    hax = gca;
    set(hax, 'Position', [pos(1), 0.05, pos(3), 0.42])
    
    set(lh, 'Box', 'off');%legend(lh, 'boxoff');
    set(lh, 'FontSize', 10);

end