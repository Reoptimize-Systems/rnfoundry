% makedocs.m
%
% Script to generate files for mbdyn docs and build the docs

% files required for embedding in docs
thisfiledir = getmfilepath ('mbdyn.makedocs');

imagedir = fullfile (thisfiledir, '..', 'images');


%% Run examples and save figures

% save_example_figs ('example_double_pendulum', imagedir, 1:6);
save_example_figs ('example_basic_cosimulation', imagedir, 1:3);


%% Make API reference
apidir = fullfile (thisfiledir, '..', 'api_reference');

mkdir (apidir);

pre_apidir = fullfile (apidir, 'pre');
mkdir (pre_apidir);

% generated docs from class file help
mbdyn_pre_class_list = dir (fullfile (getmfilepath ('mbdyn.pre.element'), '*.m'));

exclude_list = { 'baseSystem.m', ...
                 'subsystem.m' };


for ind = 1:numel(mbdyn_pre_class_list)
    
    if ~any (strcmp (mbdyn_pre_class_list(ind).name, exclude_list )) 
    
        rst = help2rst (['mbdyn.pre.', mbdyn_pre_class_list(ind).name(1:end-2)], ...
                        'SectionChars', '=-^"+');

        outpufile = fullfile (pre_apidir, [mbdyn_pre_class_list(ind).name(1:end-2), '.rst']);

        cellstr2txtfile (outpufile, rst);
    
    end
    
end


mint_apidir = fullfile (apidir, 'mint');
mkdir (mint_apidir);

% generated docs from class file help
mbdyn_pre_class_list = dir (fullfile (getmfilepath ('mbdyn.mint.MBCNodal'), '*.m'));

exclude_list = { 'externalForces.m', ...
                 'torque.m', ...
                 'force.m' };

for ind = 1:numel(mbdyn_pre_class_list)
    
    if ~any (strcmp (mbdyn_pre_class_list(ind).name, exclude_list )) 
    
        rst = help2rst (['mbdyn.mint.', mbdyn_pre_class_list(ind).name(1:end-2)], ...
                        'SectionChars', '=-^"+');

        outpufile = fullfile (mint_apidir, [mbdyn_pre_class_list(ind).name(1:end-2), '.rst']);

        cellstr2txtfile (outpufile, rst);
    
    end
    
end

post_apidir = fullfile (apidir, 'post');
mkdir (post_apidir);

rst = help2rst ( 'mbdyn.postproc', ...
                 'SectionChars', '=-^"+');

outpufile = fullfile (post_apidir, 'mbdyn.postproc.rst');

cellstr2txtfile (outpufile, rst);
        
        
%% build the docs
cd (fullfile (apidir, '..'));

[status, out] = cleansystem ('make html');

%% functions

function save_example_figs (exname, imdir, fignums_to_save)

    close all
    run (exname);
    
    handles = findall(0,'type','figure');
    
    fname = fullfile (imdir, [exname, '_%d.png']);
    
    for ind = 1:numel (handles)
        
        fignum = get (handles(ind), 'Number');
        
%         tightfig (handles(ind));
        
        switch fignum
            
            case num2cell ( fignums_to_save )
                
                saveas (handles(ind), sprintf (fname, fignum));
                
            otherwise
                
        end

    end
    
    close all
    
end
