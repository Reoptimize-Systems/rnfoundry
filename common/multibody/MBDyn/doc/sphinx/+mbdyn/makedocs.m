% makedocs.m
%
% Script to generate files for mbdyn docs

% files required for embedding in docs




thisfiledir = getmfilepath ('mbdyn.makedocs');

outputdir = fullfile (thisfiledir, '..', 'api_reference');

mkdir (outputdir);

% generated docs from class file help
mbdyn_pre_class_list = dir (fullfile (getmfilepath ('mbdyn.pre.element'), '*.m'));

exclude_list = { 'baseSystem.m' };


for ind = 1:numel(mbdyn_pre_class_list)
    
    if ~any (strcmp (mbdyn_pre_class_list(ind).name, exclude_list )) 
    
        rst = help2rst (['mbdyn.pre.', mbdyn_pre_class_list(ind).name(1:end-2)], ...
                        'SectionChars', '=-^"+');

        outpufile = fullfile (outputdir, [mbdyn_pre_class_list(ind).name(1:end-2), '.rst']);

        cellstr2txtfile (outpufile, rst);
    
    end
    
end

cd (fullfile (outputdir, '..'));

[status, out] = cleansystem ('make html');