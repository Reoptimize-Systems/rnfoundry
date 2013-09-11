function legendstrings = plotforces_linear(T, results, skip)
% plotforces_linear: plots the various forces that might arise in a linear
% machine simulation of various types

    if nargin < 3 
        skip = 1;
    end
    
    % Plot the forces and acceleration on the same graph, starting
    % with the forces
    legendstrings = {};
    currentylim = [0,0];
    
    FTtotal = zeros(size(T));
    
    if isfield(results, 'Fpto')
        plot(T(1:skip:length(T)), results.Fpto, ':b');
        currentylim = [-1.1*max(abs(results.Fpto)), 1.1*max(abs(results.Fpto))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'PTO Force'}];
        
        FTtotal = FTtotal + results.Fpto;
    end
    
    if isfield(results, 'Feff')
        plot(T(1:skip:length(T)), results.Feff, '--b');
        currentylim = [-1.1*max(abs(results.Feff)), 1.1*max(abs(results.Feff))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Effector Forces'}];
        
        FTtotal = FTtotal + results.Feff;
    end
    
    if isfield(results, 'Freac')
        plot(T(1:skip:length(T)), results.Fpto, '--b');
        currentylim = [-1.1*max(abs(results.Fpto)), 1.1*max(abs(results.Fpto))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Reactor Forces'}];
        
        FTtotal = FTtotal + results.Fpto;
    end
    
    if isfield(results, 'Ffea')
        plot(T(1:skip:length(T)), results.Ffea, '--b');
        currentylim = [-1.1*max(abs(results.Ffea)), 1.1*max(abs(results.Ffea))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'EM Forces on Translator'}];

        FTtotal = FTtotal + results.Ffea;
    end
    
    % if present plot the spring forces
    if isfield(results, 'Fs')
        hold on
        plot(T(1:skip:length(T)), results.Fs, ':g');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.Fs))), max(currentylim(2), 1.1*max(abs(results.Fs)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Spring Force'}];
    end
    
    if isfield(results, 'FfA')
        hold on
        plot(T(1:skip:length(T)), results.FfA, 'r');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FfA))), max(currentylim(2), 1.1*max(abs(results.FfA)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Armature Friction'}];
    end
    
    if isfield(results, 'FfT')
        hold on
        plot(T(1:skip:length(T)), results.FfT, 'r');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FfT))), max(currentylim(2), 1.1*max(abs(results.FfT)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Field Friction'}];
        
        FTtotal = FTtotal + results.FfT;
    end
    
    if isfield(results, 'FBuoyFric')
        hold on
        plot(T(1:skip:length(T)), results.FBuoyFric, 'r');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FBuoyFric))), max(currentylim(2), 1.1*max(abs(results.FBuoyFric)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Buoy System Friction'}];
        
        FTtotal = FTtotal + results.FBuoyFric;
    end
    
    if isfield(results, 'Fsnap')
        hold on
        plot(T(1:skip:length(T)), results.Fsnap, 'k');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.Fsnap))), max(currentylim(2), 1.1*max(abs(results.Fsnap)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Snapping Force (Reactor)'}];
        
        FTtotal = FTtotal - results.Fsnap;
        
    end
    
    if isfield(results, 'FdragF')
        hold on
        plot(T(1:skip:length(T)), results.FdragF, ':r');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FdragF))), max(currentylim(2), 1.1*max(abs(results.FdragF)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Fluid Drag (Reactor)'}];
        
%         FTtotal = FTtotal - results.FdragF;
        
    end
    
    if isfield(results, 'FdragA')
        hold on
        plot(T(1:skip:length(T)), results.FdragA, ':r');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FdragA))), max(currentylim(2), 1.1*max(abs(results.FdragA)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Fluid Drag (Reactor)'}];
        
%         FTtotal = FTtotal - results.FdragA;
        
    end
    
    if isfield(results, 'FLinearDrag')
        hold on
        plot(T(1:skip:length(T)), results.FLinearDrag, 'm');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FLinearDrag))), max(currentylim(2), 1.1*max(abs(results.FLinearDrag)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Linear Drag Force (Reactor)'}];
        
        FTtotal = FTtotal - results.FLinearDrag;
    end
    
    if isfield(results, 'FLoss')
        hold on
        plot(T(1:skip:length(T)), results.FLoss, '--m');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FLoss))), max(currentylim(2), 1.1*max(abs(results.FLoss)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Loss Force (Effector)'}];
        
        FTtotal = FTtotal - results.FLoss;
    end
    
    if isfield(results, 'Fes')
        hold on
        plot(T(1:skip:length(T)), results.Fes, 'g');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.Fes))), max(currentylim(2), 1.1*max(abs(results.Fes)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'End Stop Force (Effector)'}];
        
        FTtotal = FTtotal - results.Fes;
    end
    
    if isfield(results, 'Faddtrans')
        
        hold on
        plot(T(1:skip:length(T)), results.Faddtrans(:,1), ':k');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.Faddtrans(:,1)))), max(currentylim(2), 1.1*max(abs(results.Faddtrans(:,1))))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Total Translator Force'}];
    end

    xlabel('Time (s)');%,'FontSize',20);
    ylabel('Force (N)');%, 'FontSize', 20);

end