function savesplitoderesult(simoptions, basename, results, sol, odeinternals)
% saves intermediate results from a split ODE sim to disk
%
% Syntax
%
% savesplitoderesult(simoptions, basename, results, sol, odeinternals)
%
% 

    if isfield(simoptions, 'SplitSaveDir')
        splitsavedir = simoptions.SplitSaveDir;
    else
        splitsavedir = '';
    end
    if results.block == 1
        s = warning('off', 'MATLAB:DELETE:FileNotFound');
        delete(fullfile(splitsavedir,[basename, '_*.mat']));
        %  warning('on', 'MATLAB:DELETE:FileNotFound');
        warning(s);
    end
    save(fullfile(splitsavedir,[basename, '_', num2str(results.block), '.mat']), 'sol', 'odeinternals');

end