function [docrootdir, zipfilename] = makedocs (varargin)
% makedocs.m
%
% Script to generate files for wave docs and build the docs

    options.Version = 'devel';
    options.ToolboxName = 'EWST';
%     options.BuildInTempDir = false;
    
    options = parse_pv_pairs (options, varargin);
    
    origpwd = pwd ();
    
    CC = onCleanup (@() cd (origpwd));

    % files required for embedding in docs
    thisfiledir = getmfilepath ('wsim.makedocs');
    
    docrootdir = fullfile (thisfiledir, '..');

%     if options.BuildInTempDir
%         
%        
%         build_docrootdir = fullfile (tempdir (), outfilename);
%         
%     else
%         
%     end
    
    %% Run examples and save figures
    
    % imagedir = fullfile (docrootdir, 'images');
    %
    % save_example_figs ('example_double_pendulum', imagedir, 1:6);
    % save_example_figs ('example_basic_cosimulation', imagedir, 1:3);


    %% Make API reference
    wsim_apidir = fullfile (docrootdir, 'api_reference');

    mkdir (wsim_apidir);

    % generated docs from class file help
    wsim_class_list = dir (fullfile (getmfilepath ('wsim.wecSim'), '*.m'));

    % exclude_list = { 'baseSystem.m', ...
    %                  'subsystem.m' };
    exclude_list = { 'mooring.m' };

    for ind = 1:numel(wsim_class_list)

        if ~any (strcmp (wsim_class_list(ind).name, exclude_list )) 

            rst = help2rst (['wsim.', wsim_class_list(ind).name(1:end-2)], ...
                             'SectionChars', '=-^"+');

            outputfile = fullfile (wsim_apidir, [wsim_class_list(ind).name(1:end-2), '.rst']);

            cellstr2txtfile (outputfile, rst);

        end

    end


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
    rmdir (temp_docs_dir);
    
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


% function save_example_figs (exname, imdir, fignums_to_save)
% 
%     close all
%     run (exname);
%     
%     handles = findall(0,'type','figure');
%     
%     fname = fullfile (imdir, [exname, '_%d.png']);
%     
%     for ind = 1:numel (handles)
%         
%         fignum = get (handles(ind), 'Number');
%         
% %         tightfig (handles(ind));
%         
%         switch fignum
%             
%             case num2cell ( fignums_to_save )
%                 
%                 saveas (handles(ind), sprintf (fname, fignum));
%                 
%             otherwise
%                 
%         end
% 
%     end
%     
%     close all
%     
% end