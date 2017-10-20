% generate_help_files.m
%

% list of help commands to output to files
helplist = { ...
'nemoh.simulation/writeNemoh', ...
'nemoh.simulation/run' ...
};


%%
thisfiledir = getmfilepath ('generate_help_files.m');

outdir = fullfile (thisfiledir, 'generated');

mkdir (outdir);


for ind = 1:numel (helplist)
    
    helptxt = help (helplist{ind});
    
    fname = fullfile (outdir, sprintf ('help.%s', strrep (helplist{ind}, '/', '.')));
    
    fid = fopen (fname, 'w');
    CC = onCleanup (@() fclose (fid));
    
    fprintf (fid, '%s', helptxt);
    
    clear CC
    
end



