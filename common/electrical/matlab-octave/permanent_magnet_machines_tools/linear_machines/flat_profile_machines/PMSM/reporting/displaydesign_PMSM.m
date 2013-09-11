function displaydesign_PMSM(design, simoptions)
% displaydesign_PMSM: displays information about the PMSM machine design
%

    fprintf(1, 'PMSM Design:\n');
    fprintf(1, 'Dim        Val      Dim        Val      Dim        Val      Dim        Val\n');
    fprintf(1, 'Wp:     %8.4f    Wm:     %8.4f    Wc:     %8.4f    Wt:     %8.4f\n', design.Wp, design.Wm, design.Wc, design.Wt);
    fprintf(1, 'Ws:     %8.4f    hba:    %8.4f    ht:     %8.4f    g:      %8.4f\n', design.Ws, design.hba, design.ht, design.g);
    fprintf(1, 'hm:     %8.4f    hbf:    %8.4f    ls:     %8.4f    FF:     %8.4f\n', design.hm, design.hbf, design.ls, design.fillfactor);

    if isfield(design, 'poles')
        fprintf(1, 'poles: %s', num2str(design.poles));
    end
    
    if isfield(design, 'Rs1')
        fprintf(1, '\tHmag: %f\t', design.Hmag);
    end       

    fprintf(1, '\n');
    
    displaydesign_AM(design, simoptions);
    
%     lastchar = '';
%     
%     if isfield(design, 'Rs1')
%         fprintf(1, 'Rs1: %f\t', design.Rs1);
%         lastchar = '\n';
%     end
%     
%     if isfield(design, 'Rs2')
%         fprintf(1, 'Rs2: %f\t', design.Rs1);
%         lastchar = '\n';
%     end
%     
%     if isfield(design, 'Ws1')
%         fprintf(1, 'Ws1: %f\t', design.Rs1);
%         lastchar = '\n';
%     end
%     
%     if isfield(design, 'Ws2')
%         fprintf(1, 'Ws2: %f\t', design.Rs1);
%         lastchar = '\n';
%     end
    
%     fprintf(1, lastchar);
    
    
end