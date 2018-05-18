function path = find_mbdyn (throw)
% locate the mbdyn executeable
%
% Syntax
%
% path = find_mbdyn ()
%
% Description
%
% findmbdyn searches in set of predetermined locations for an MBDyn
% executable. It first searches the current directory, then various
% directories in which it might be expected to be found on different
% systems.
%
% Input
%
%  throw - flag determining whether to throw an error if mbdyn is not found
%
% Output
%
%  path - path to the MBDyn executeable which was found
%
%
% See Also: mbdyn.mint.start_mbdyn, mbdyn.mint.find_libmbc
%

    if nargin < 1
        throw = true;
    end
    
    mbdyn.pre.base.checkLogicalScalar (throw, true, 'throw');

    candidate_locs = { ...
        fullfile(pwd (), 'mbdyn'), ...
        fullfile(pwd (), 'mbdyn.exe'), ...
        'c:\Program Files (x86)\MBDyn\mbdyn.exe', ...
        'c:\Program Files\MBDyn\mbdyn.exe', ...
        '/usr/local/mbdyn/bin/mbdyn', ...
        '/usr/local/bin/mbdyn' ...
        '/opt/mbdyn/bin/mbdyn', ...
        '/opt/bin/mbdyn' ...
        '/opt/mbdyn' ...
                     };

    comparch = computer ('arch');
    switch comparch

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
        switch comparch

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

    path = '';
    for ind = 1:numel (candidate_locs)
        if exist (candidate_locs{ind}, 'file') == 2
            path = candidate_locs{ind};
            return;
        end
    end

    if throw
        error ('The MBDyn executeable could not be found.');
    end

end