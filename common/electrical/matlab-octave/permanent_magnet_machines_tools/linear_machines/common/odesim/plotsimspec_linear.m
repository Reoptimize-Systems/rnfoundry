function plotsimspec_linear (simoptions)

    figure; 
    
    [ hAx, hp1, hp2 ] = plotyy (simoptions.drivetimes, simoptions.xT, ...
            simoptions.drivetimes, simoptions.vT );
%         
%     
    
    
    hold (hAx(1), 'on');
    hold (hAx(2), 'on');
    
    t = linspace (simoptions.drivetimes(1), simoptions.drivetimes(end), 10*numel(simoptions.drivetimes));
    
    [xTinterp, vTinterp] = prescribedmotvelpos (t, simoptions);
    
    hp1pp = plot (hAx(1), t, xTinterp, ':');
    hp2pp = plot (hAx(2), t, vTinterp, ':');
    
    
    hp1.LineStyle = 'none';
    hp1.Marker = 'x';
    p1col = hp1.Color;
    
    hp2.LineStyle = 'none';
    hp2.Marker = 'x';
    p2col = hp2.Color;
    
    hp1pp.Color = p1col;
    hp2pp.Color = p2col;
    
    hold off;
    
    legend ('xT', 'xT Interp' ,'vT' ,'vT Interp', 'Location', 'Best');
    
    xlabel ('Time [s]');
    
    ylabel (hAx(1), 'Position [m]');
    ylabel (hAx(2), 'Velocity [ms^{-1}]');
    
    hold (hAx(1), 'off');
    hold (hAx(2), 'off');
    
end
