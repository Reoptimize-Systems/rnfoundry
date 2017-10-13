function path = findmbdyn ()
% locate the mbdyn executeable

    candidate_locs = { ...
        fullfile(pwd (), 'mbdyn'), ...
        fullfile(pwd (), 'mbdyn.exe'), ...
        'c:\Program Files (x86)\MBDyn\mbdyn.exe', ...
        'c:\Program Files\MBDyn\mbdyn.exe', ...
        '/usr/local/bin/mbdyn/mbdyn', ...
        '/usr/local/bin/mbdyn' ...
                     };

    switch computer ('arch')

        case 'win64'
            candidate_locs = [candidate_locs, ...
                fullfile(getmfilepath ('mexmbdyn_setup'), 'x86_64-w64-mingw32', 'bin', 'mbdyn.exe')];
        case 'win32'
            candidate_locs = [candidate_locs, ...
                fullfile(getmfilepath ('mexmbdyn_setup'), 'i686-w64-mingw32', 'bin', 'mbdyn.exe')];
        case 'glnxa64'
            candidate_locs = [candidate_locs, ...
                fullfile(getmfilepath ('mexmbdyn_setup'), 'x86_64-linux-gnu', 'bin', 'mbdyn')];
        case 'glnxa32'
            candidate_locs = [candidate_locs, ...
                fullfile(getmfilepath ('mexmbdyn_setup'), 'i686-linux-gnu', 'bin', 'mbdyn')];

    end

    if exist ('rnfoundry_setup', 'file') == 2
        switch computer ('arch')

            case 'win64'
                candidate_locs = [candidate_locs, ...
                    fullfile(getmfilepath ('rnfoundry_setup'), 'x86_64-w64-mingw32', 'bin', 'mbdyn.exe')];
            case 'win32'
                candidate_locs = [candidate_locs, ...
                    fullfile(getmfilepath ('rnfoundry_setup'), 'i686-w64-mingw32', 'bin', 'mbdyn.exe')];
            case 'glnxa64'
                candidate_locs = [candidate_locs, ...
                    fullfile(getmfilepath ('rnfoundry_setup'), 'x86_64-linux-gnu', 'bin', 'mbdyn')];
            case 'glnxa32'
                candidate_locs = [candidate_locs, ...
                    fullfile(getmfilepath ('rnfoundry_setup'), 'i686-linux-gnu', 'bin', 'mbdyn')];

        end
    end

    for ind = 1:numel (candidate_locs)
        if exist (candidate_locs{ind}, 'file') == 2
            path = candidate_locs{ind};
            return;
        end
    end

    error ('The MBDyn executeable could not be found.');

end