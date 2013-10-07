function displaydesign_ACPMSM(design, simoptions)
% displays information about the ACPMSM machine design
%


    fprintf(1, 'ACPMSM Design:\n');
    fprintf(1, 'Dim        Val      Dim        Val      Dim        Val      Dim        Val\n');
    fprintf(1, 'lm:     %8.4f    Taup:   %8.4f    bp:     %8.4f    dg:     %8.4f\n', design.lm, design.Taup, design.bp, design.dg);
    fprintf(1, 'ls:     %8.4f    dbi:    %8.4f    Wc:     %8.4f    Hc:     %8.4f\n', design.ls, design.dbi, design.Wc, design.Hc);
    fprintf(1, 'g:      %8.4f\n', design.g, design.g);

    if isfield(design, 'Poles')
        fprintf(1, 'Poles: %s\n', num2str(design.Poles));
    end

    displaydesign_AM(design, simoptions)
    
end