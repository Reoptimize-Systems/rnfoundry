% pkg_gapolyfitn

gapolyfitndir = fileparts(which('gapolyfitn'));

polyfitndir = fileparts(which('polyfitn'));

file_list = getFileDependencies(('gapolyfitn_example_script.m'));

expdir = fullfile(gapolyfitndir, 'gapolyfitndir');

if ~exist(expdir, 'dir')
    mkdir(expdir);
end

if ~exist(fullfile(expdir, 'private'), 'dir')
    mkdir(expdir, 'private')
end

for i = 1:numel(file_list)

    if ~isempty(findstr(file_list{i}, fullfile('matlab', 'GA_Toolbox', 'genetic'))) || ...
            ~isempty(findstr(file_list{i}, fullfile('matlab', 'Multicore')))
        
        % leave it
        fprintf(1, 'ignoring: %s\n', file_list{i})
        
    elseif ~isempty(findstr(file_list{i}, fullfile('matlab', 'PolyfitnTools', 'gapolyfitn')))
        
        fprintf(1, 'copying: %s to %s\n', file_list{i}, expdir)
        
        copyfile(file_list{i}, expdir);
        
    else
        
        fprintf(1, 'copying: %s to %s\n', file_list{i}, fullfile(expdir, 'private'))
        
        copyfile(file_list{i}, fullfile(expdir, 'private'));
        
    end

end

zip(fullfile(gapolyfitndir, 'gapolyfitn.zip'), expdir);

rmdir(expdir, 's')





