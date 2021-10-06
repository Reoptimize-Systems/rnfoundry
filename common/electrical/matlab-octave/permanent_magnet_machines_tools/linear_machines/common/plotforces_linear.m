function legendstrings = plotforces_linear (T, results, skip)
% plots the various forces that might arise in a linear machine ODE
% simulation of various types
%
% Syntax
%
% legendstrings = plotforces_linear (T, results, skip)
%
% Input
%
%  T - time series as output by the ode solver
%
%  results - structure containing one or more of the following fields:
%     
%    Fpto
%    Feff
%    Freac
%    Ffea
%    Fs
%    FfA
%    FfT
%    FBuoyFric
%    FdragF
%    FdragA
%    FLinearDrag
%    FLoss
%    Fes
%    Faddtrans
%    Fsnap
%
%  skip - optional scalar determining how much output to produce, indices
%    1:skip:end of the solution vectors will be plotted
%
%

    if nargin < 3 
        skip = 1;
    end
    
    % Plot the forces and acceleration on the same graph, starting
    % with the forces
    legendstrings = {};
    currentylim = [0,0];
    
    FTtotal = zeros( size ( T(1:skip:length(T)) ) );
    
    if isfield(results, 'Fpto')
        plot(T(1:skip:length(T)), results.Fpto(1:skip:length(T),:), ':b');
        currentylim = [-1.1*max(abs(results.Fpto(1:skip:length(T),:))), 1.1*max(abs(results.Fpto(1:skip:length(T),:)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'PTO Force'}];
        
        FTtotal = FTtotal + results.Fpto;
    end
    
    if isfield(results, 'Feff')
        plot(T(1:skip:length(T)), results.Feff(1:skip:length(T),:), '--b');
        currentylim = [-1.1*max(abs(results.Feff(1:skip:length(T),:))), 1.1*max(abs(results.Feff(1:skip:length(T),:)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Effector Forces'}];
        
        FTtotal = FTtotal + results.Feff(1:skip:length(T),:);
    end
    
    if isfield(results, 'Freac')
        plot(T(1:skip:length(T)), results.Freac(1:skip:length(T),:), '--b');
        currentylim = [-1.1*max(abs(results.Freac(1:skip:length(T),:))), 1.1*max(abs(results.Freac(1:skip:length(T),:)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Reactor Forces'}];
        
        FTtotal = FTtotal + results.Freac(1:skip:length(T),:);
    end
    
    if isfield(results, 'Ffea')
        plot(T(1:skip:length(T)), results.Ffea(1:skip:length(T),:), '--b');
        currentylim = [-1.1*max(abs(results.Ffea(1:skip:length(T),:))), 1.1*max(abs(results.Ffea(1:skip:length(T),:)))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'EM Forces on Translator'}];

        FTtotal = FTtotal + results.Ffea(1:skip:length(T),:);
    end
    
    % if present plot the spring forces
    if isfield(results, 'Fs')
        hold on
        plot(T(1:skip:length(T)), results.Fs(1:skip:length(T),:), ':g');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.Fs(1:skip:length(T),:)))), max(currentylim(2), 1.1*max(abs(results.Fs(1:skip:length(T),:))))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Spring Force'}];
    end
    
    if isfield(results, 'FfA')
        hold on
        plot(T(1:skip:length(T)), results.FfA(1:skip:length(T),:), 'r');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FfA(1:skip:length(T),:)))), max(currentylim(2), 1.1*max(abs(results.FfA(1:skip:length(T),:))))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Armature Friction'}];
    end
    
    if isfield(results, 'FfT')
        hold on
        plot(T(1:skip:length(T)), results.FfT, 'r');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FfT(1:skip:length(T),:)))), max(currentylim(2), 1.1*max(abs(results.FfT(1:skip:length(T),:))))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Field Friction'}];
        
        FTtotal = FTtotal + results.FfT(1:skip:length(T),:);
    end
    
    if isfield(results, 'FBuoyFric')
        hold on
        plot(T(1:skip:length(T)), results.FBuoyFric(1:skip:length(T),:), 'r');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FBuoyFric(1:skip:length(T),:)))), max(currentylim(2), 1.1*max(abs(results.FBuoyFric(1:skip:length(T),:))))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Buoy System Friction'}];
        
        FTtotal = FTtotal + results.FBuoyFric(1:skip:length(T),:);
    end
    
    if isfield(results, 'Fsnap')
        hold on
        plot(T(1:skip:length(T)), results.Fsnap(1:skip:length(T),:), 'k');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.Fsnap(1:skip:length(T),:)))), max(currentylim(2), 1.1*max(abs(results.Fsnap(1:skip:length(T),:))))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Snapping Force (Reactor)'}];
        
        FTtotal = FTtotal - results.Fsnap(1:skip:length(T),:);
        
    end
    
    if isfield(results, 'FdragF')
        hold on
        plot(T(1:skip:length(T)), results.FdragF(1:skip:length(T),:), ':r');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FdragF(1:skip:length(T),:)))), max(currentylim(2), 1.1*max(abs(results.FdragF(1:skip:length(T),:))))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Fluid Drag (Reactor)'}];
        
%         FTtotal = FTtotal - results.FdragF;
        
    end
    
    if isfield(results, 'FdragA')
        hold on
        plot(T(1:skip:length(T)), results.FdragA(1:skip:length(T),:), ':r');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FdragA(1:skip:length(T),:)))), max(currentylim(2), 1.1*max(abs(results.FdragA(1:skip:length(T),:))))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Fluid Drag (Reactor)'}];
        
%         FTtotal = FTtotal - results.FdragA;
        
    end
    
    if isfield(results, 'FLinearDrag')
        hold on
        plot(T(1:skip:length(T)), results.FLinearDrag(1:skip:length(T),:), 'm');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FLinearDrag(1:skip:length(T),:)))), max(currentylim(2), 1.1*max(abs(results.FLinearDrag(1:skip:length(T),:))))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Linear Drag Force (Reactor)'}];
        
        FTtotal = FTtotal - results.FLinearDrag(1:skip:length(T),:);
    end
    
    if isfield(results, 'FLoss')
        hold on
        plot(T(1:skip:length(T)), results.FLoss(1:skip:length(T),:), '--m');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.FLoss(1:skip:length(T),:)))), max(currentylim(2), 1.1*max(abs(results.FLoss(1:skip:length(T),:))))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Loss Force (Effector)'}];
        
        FTtotal = FTtotal - results.FLoss(1:skip:length(T),:);
    end
    
    if isfield(results, 'Fes')
        hold on
        plot(T(1:skip:length(T)), results.Fes(1:skip:length(T),:), 'g');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.Fes(1:skip:length(T),:)))), max(currentylim(2), 1.1*max(abs(results.Fes(1:skip:length(T),:))))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'End Stop Force (Effector)'}];
        
        FTtotal = FTtotal - results.Fes(1:skip:length(T),:);
    end
    
    if isfield(results, 'Faddtrans')
        
        hold on
        plot(T(1:skip:length(T)), results.Faddtrans(:,1), ':k');
        hold off
        currentylim = [min(currentylim(1), -1.1*max(abs(results.Faddtrans(1:skip:length(T),1)))), max(currentylim(2), 1.1*max(abs(results.Faddtrans(1:skip:length(T),1))))];
        ylim(currentylim);
        legendstrings = [legendstrings; {'Total Translator Force'}];
    end

    xlabel('Time (s)');%,'FontSize',20);
    ylabel('Force (N)');%, 'FontSize', 20);

end