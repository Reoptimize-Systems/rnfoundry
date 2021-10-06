function displaydesign_ACTIAM(design, simoptions)
% displaydesign_ACTM: displays information about the ACTM machine design
%

    fprintf(1, 'ACTIAM Design:\n');
    fprintf(1, 'Dim        Val      Dim        Val      Dim        Val      Dim        Val\n');
    fprintf(1, 'Rm:     %8.4f    Ri:     %8.4f    Ro:     %8.4f    Ra:     %8.4f\n', design.Rm, design.Ri, design.Ro, design.Ra);
    fprintf(1, 'Rso:    %8.4f    Rsi:    %8.4f    Wp:     %8.4f    Wm:     %8.4f\n', design.Rso, design.Rsi, design.Wp, design.Wm);
    fprintf(1, 'Ws:     %8.4f    Wc      %8.4f    Hc      %8.4f    g:      %8.4f\n', design.Ws, design.Wc, design.Hc, design.g);
    fprintf(1, 'Hmag:   %8.4f\n', design.Hmag);

    if isfield(design, 'Poles')
        fprintf(1, 'Poles: %s', num2str(design.Poles));
    end

    fprintf(1, '\n');
    
    lastchar = '';
    
    if isfield(design, 'Rs1')
        fprintf(1, 'Rs1: %f\t', design.Rs1);
        lastchar = '\n';
    end
    
    if isfield(design, 'Rs2')
        fprintf(1, 'Rs2: %f\t', design.Rs1);
        lastchar = '\n';
    end
    
    if isfield(design, 'Ws1')
        fprintf(1, 'Ws1: %f\t', design.Rs1);
        lastchar = '\n';
    end
    
    if isfield(design, 'Ws2')
        fprintf(1, 'Ws2: %f\t', design.Rs1);
        lastchar = '\n';
    end
    
    fprintf(1, lastchar);
    
    displaydesign_AM(design, simoptions)
    
    
end