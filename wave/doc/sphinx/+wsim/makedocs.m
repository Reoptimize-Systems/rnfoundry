% makedocs.m
%
% Script to generate files for mbdyn docs and build the docs

TNShort = 'EWST';

% files required for embedding in docs
thisfiledir = getmfilepath ('wsim.makedocs');

imagedir = fullfile (thisfiledir, '..', 'images');


%% Run examples and save figures

% save_example_figs ('example_double_pendulum', imagedir, 1:6);
% save_example_figs ('example_basic_cosimulation', imagedir, 1:3);


%% Make API reference
wsim_apidir = fullfile (thisfiledir, '..', 'api_reference');

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
cd (fullfile (wsim_apidir, '..'));

[status, out] = cleansystem ('make html');

%% export docs to a zip file

% [status, out] = cleansystem ('make html');

%% functions
% 
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
