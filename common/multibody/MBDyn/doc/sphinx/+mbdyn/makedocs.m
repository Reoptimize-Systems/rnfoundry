function [ docrootdir, zipfilename ] = makedocs (varargin)
% makedocs.m
%
% function to generate files for mbdyn docs and build the docs

    options.Version = 'devel';
    options.ToolboxName = 'mbdyn-matlab-toolbox';
%     options.BuildInTempDir = false;
    
    options = parse_pv_pairs (options, varargin);
    
    origpwd = pwd ();
    
    CC = onCleanup (@() cd (origpwd));% makedocs.m

    % files required for embedding in docs
    thisfiledir = getmfilepath ('mbdyn.makedocs');
    
    docrootdir = fullfile (thisfiledir, '..');

    imagedir = fullfile (docrootdir, 'images');

    %% Run examples and save figures

    % save_example_figs ('example_double_pendulum', imagedir, 1:6);
    save_example_figs ('example_basic_cosimulation', imagedir, 1:3);

    %% Make API reference
    apidir = fullfile (docrootdir, 'api_reference');

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
    cd (fullfile (docrootdir));

    [status, out] = cleansystem ('make html');

    disp (out);
    
    if status ~= 0
        error ('Doc build failed');
    end
    
    %% export docs to a zip file
    
    html_dir = fullfile (docrootdir, '_build', 'html');

    outfilename = sprintf ('%s-%s-html-docs', options.ToolboxName, options.Version);
    
    % create a directory name for the docs
    temp_docs_dir = fullfile (tempdir (), outfilename);
    
    % make sure this doesn't exist already
    delete (temp_docs_dir);
    
    copyfile (html_dir, temp_docs_dir);
    
    cd (tempdir ());
    
    [status, out] = cleansystem (sprintf ('zip -qr %s.zip %s/', outfilename, outfilename));
    
    disp (out);
    
    zipfilename = [outfilename, '.zip'];
    
    if status == 0
        copyfile (fullfile (tempdir (), zipfilename), fullfile (docrootdir, '..'));
    else
        error ('Making zip failed');
    end
    
    
end


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
