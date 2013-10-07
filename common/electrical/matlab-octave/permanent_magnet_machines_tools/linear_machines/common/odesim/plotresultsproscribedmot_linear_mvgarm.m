function plotresultsproscribedmot_linear_mvgarm(T, Y, results, design, skip)
% plots the output of a simulation performed using the
% proscribedmot_linear_mvgarm function
%
% Syntax
% 
% plotresultsproscribedmot_linear_mvgarm(T, Y, results, design)
% plotresultsproscribedmot_linear_mvgarm(..., skip)
%
% Inputs
%
%  T - vector of time steps for which the values of the the snapper
%      outputs in Y have been evaluated
%
%  Y - This is the output from the ode solver for a proscribed motion
%      simulation. In this case the vector contains the armature
%      displacment, velocity and coil current respectively at each time
%      step in the vector T. So that:
%
%       Y = [ Iphases, xA, vA ]
%
%  results - A structure containing some further outputs from the
%            simulation. It should contain the following fields: 
%
%            EMF: a vector of the EMF values at each time step in the
%             vector T. 
%
%            Ffea: the electromagnetic forces at each time step in T
%
%            Fs: the spring forces at each time step in T
%
%            ydot: an (n x p) matrix containing the results of the 
%             right hand side of the differential equations solved by
%             the matlab ode solvers.
%
%            The following fields can also optionally be supplied:
%
%            FDrive: The drive force applied to the translator at each
%             time step in T.
%
%            xT: Translator position
%
%            vT: Translator velocity
%
%            Ff: (optional) The frictional forces between the armature and
%             translator at each time step in T.
%
%   skip - optional value determining what proportion of time points are to be
%          plotted. If skip is 2, every other step is plotted, if 3, every
%          third point etc. Defaults to 1 if not supplied (no solution
%          points skipped)

    % Y is the coil current for each phase, where each phase has it's own
    % column
    
    if nargin < 5
        skip = 1;
    end
    
    Icoils = Y(:,1:design.Phases);
    
    xA = Y(:,design.Phases + 1);
    vA = Y(:,design.Phases + 2);
    
%     
%     % Get the Y-axis limits for the various plots
%     maxY = [1.05*max([abs([Icoils, results.xT, results.vT, xA, vA]); repmat(1e-10,1,size(Icoils,2)+4)], [], 1)];

    % create the figure, ensuring it is visible
    h = figure('visible','on', 'Units','normalized','outerposition',[0 0.02 1 0.98]);
    
    % we will create a 2 x 1 subplot
    subplot(2,1,1);
    
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
                                                     'xA', xA, ...
                                                     'vA', vA, ...
                                                     'AxesPosition', [0.14/2, 0.55, 0.86, 0.42]);
                                                    
    lh = legend(legendstrings, 'Location','NorthWest');

    set(lh, 'Box', 'off');
    set(lh, 'FontSize', 10);

    pos = get(hax, 'Position');

    % create a subplot for the forces in the results structure
    subplot(2,1,2);

    legendstrings = plotforces_linear(T, results, skip);
    
    lh = legend(legendstrings, 'Location', 'NorthWest');
    
    hax = gca;
    set(hax, 'Position', [pos(1), 0.05, pos(3), 0.42])
    
    set(lh, 'Box', 'off');%legend(lh, 'boxoff');
    set(lh, 'FontSize', 10);

end