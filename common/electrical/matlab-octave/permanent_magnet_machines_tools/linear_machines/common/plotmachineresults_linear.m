function [hax, legendstrings] = plotmachineresults_linear(T, varargin)

    Inputs.xF = [];
    Inputs.vF = [];
    Inputs.xA = [];
    Inputs.vA = [];
    Inputs.Iphases = [];
    Inputs.EMF = [];
    Inputs.xT = [];
    Inputs.vT = [];
    Inputs.AxesPosition = [0.14/2, 0.55, 0.86, 0.42];
    Inputs.hax = [];
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isempty(Inputs.hax)
        hax = axes;
    else
        hax = Inputs.hax;
    end
    
    set(hax, 'Position', Inputs.AxesPosition)
    
    legendstrings = {};
    
    if ~isempty(Inputs.xA)
        
        % First plot the positions on the top subplot
        plot(T(:), Inputs.xA, ':b', T(:), Inputs.xT, ':r', T(:), Inputs.xT - Inputs.xA, ':m');
        ylim([-max([max(abs(Inputs.xA)), max(abs(Inputs.xT)), max(abs(Inputs.xT - Inputs.xA))]), ...
               max([max(abs(Inputs.xA)), max(abs(Inputs.xT)), max(abs(Inputs.xT - Inputs.xA))])]);
        ylabel('Displacement, m');%,'FontSize',20);
        xlabel('Time (s)');%,'FontSize',20);
        
        legendstrings = [legendstrings, {'xA : Armature Position, m', 'xT : Translator Position, m', 'xR : Relative Position, m'}];
        
    elseif ~isempty(Inputs.xF)
        
        % First plot the positions on the top subplot
        plot(T(:), Inputs.xF, ':b', T(:), Inputs.xT, ':r', T(:), Inputs.xF - Inputs.xT, ':m');
        
        ylim([-max([max(abs(Inputs.xF)), max(abs(Inputs.xT)), max(abs(Inputs.xF - Inputs.xT))]), ...
               max([max(abs(Inputs.xF)), max(abs(Inputs.xT)), max(abs(Inputs.xF - Inputs.xT))])]);
        ylabel('Displacement, m');%,'FontSize',20);
        xlabel('Time (s)');%,'FontSize',20);
        
        legendstrings = [legendstrings, {'xF : Field Position, m', 'xT : Translator Position, m', 'xR : Relative Position, m'}];
        
    else
        
        % First plot the positions on the top subplot
        plot(T(:), Inputs.xT, '--r');
        
        ylim([-max(abs(Inputs.xT)), max(abs(Inputs.xT))]);
        ylabel('Displacement, m');%,'FontSize',20);
        xlabel('Time (s)');%,'FontSize',20);
        
        legendstrings = [legendstrings, {'xT : Translator Position, m'}];
        
    end

    if ~isempty(Inputs.vA)
        
        % Now plot the velocities on their own axes
        addaxis(T(:), ...
                Inputs.vA, ...
                [-max([max(abs(Inputs.vA)), max(abs(Inputs.vT)), max(abs(Inputs.vT - Inputs.vA))]), ...
                  max([max(abs(Inputs.vA)), max(abs(Inputs.vT)), max(abs(Inputs.vT - Inputs.vA))])], ...
                '-b');
        addaxisplot(T(:), Inputs.vT, 2, '-r');
        addaxisplot(T(:), Inputs.vT - Inputs.vA, 2, '-m');
        addaxislabel(2,'Velocity, ms^{-1}');

        legendstrings = [legendstrings, {'vA : Armature Velocity, ms^{-1}', 'vT : Translator Velocity, ms^{-1}', 'vR : Relative Velocity, ms^{-1}'}];
        
    elseif ~isempty(Inputs.xF)
        
        addaxis(T(:), ...
                Inputs.vF, ...
                [-max([max(abs(Inputs.vF)), max(abs(Inputs.vT)), max(abs(Inputs.vF - Inputs.vT))]), ...
                  max([max(abs(Inputs.vF)), max(abs(Inputs.vT)), max(abs(Inputs.vF - Inputs.vT))])], ...
                '-b');
        addaxisplot(T(:), Inputs.vT, 2, '-r');
        addaxisplot(T(:), Inputs.vF - Inputs.vT, 2, '-m');
        addaxislabel(2,'Velocity, ms^{-1}');
     
        legendstrings = [legendstrings, {'vF : Field Velocity, ms^{-1}', 'vT : Translator Velocity, ms^{-1}', 'vR : Relative Velocity, ms^{-1}'}];
        
    else
        
        addaxis(T(:), ...
                Inputs.vT, ...
                [-max(abs(Inputs.vT)), max(abs(Inputs.vT))], ...
                '-r');
            
        addaxislabel(2,'Velocity, ms^{-1}');
        
        legendstrings = [legendstrings, {'vT : Translator Velocity, ms^{-1}'}];
        
    end
    
    % Now plot the current on its own axes
    addaxis(T(:), Inputs.Iphases(:,1), [-max(abs([Inputs.Iphases(:,1); 1e-10])), max(abs([Inputs.Iphases(:,1); 1e-10]))], '-k');
    addaxislabel(3, 'Current, A');

    legendstrings = [legendstrings, {'I : Phase Current, A'}];
    
    % Now plot the EMF on its own axes
    addaxis(T(:), Inputs.EMF(:,1),  [-max(abs([Inputs.EMF(:,1); 1e-10])), max(abs([Inputs.EMF(:,1); 1e-10]))], '-','color', rgb('DarkGreen'));
    addaxislabel(4,'EMF, V')

    legendstrings = [legendstrings, {'EMF, V'}];

%     set(gca, 'Position', [0.14/2, 0.55, 0.86, 0.42])

end