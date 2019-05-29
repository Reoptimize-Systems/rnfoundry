function [libdir, includedir, libwasfound, headerwasfound] = find_libmbc (varargin)
% locate the mbdyn mbc library directory and header directory
%
% Syntax
%
% [libdir, includedir, libwasfound, headerwasfound] = find_libmbc ()
%
% Description
%
% find_libmbc searches in set of predetermined locations for the MBDyn mbc
% library and header files.
%
% Output
%
%  libdir - path to the directoy where the libmbc.so or libmbc.a library
%   was found
%
%  includedir - path to the directoy where the mbc.h library header file
%   was found
%
%  libwasfound - true/false flag indicating whether the library was found
%   in any of the searched locations
%
%  headerwasfound - true/false flag indicating whether the header mbc files
%   were found in any of the searched locations
%
%
% See Also: mbdyn.mint.start_mbdyn, mbdyn.mint.find_mbdyn
%

    options.Throw = true;
    options.ForceLocalDirsOnlyForArch = '';
    
    options = parse_pv_pairs (options, varargin);
    
    mbdyn.pre.base.checkLogicalScalar (options.Throw, true, 'throw');
    if ~isempty (options.ForceLocalDirsOnlyForArch)
        mbdyn.pre.base.checkAllowedStringInputs (options.ForceLocalDirsOnlyForArch, {'win64', 'win32', 'glnxa64', 'glnxa32'} , true, 'ForceLocalDirsOnlyForArch');
    end
    
    if isempty (options.ForceLocalDirsOnlyForArch)

        libdir_candidate_locs = { ...
            'c:\Program Files (x86)\MBDyn\', ...
            'c:\Program Files (x86)\MBDyn\lib', ...
            'c:\Program Files\MBDyn\', ...
            'c:\Program Files\MBDyn\lib', ...
            '/usr/local/lib' ...
            '/usr/local/mbdyn/lib', ...
            '/opt/lib', ...
            '/opt/mbdyn/lib' ...
                         };

        includedir_candidate_locs = { ...
            'c:\Program Files (x86)\MBDyn\', ...
            'c:\Program Files (x86)\MBDyn\include', ...
            'c:\Program Files\MBDyn\', ...
            'c:\Program Files\MBDyn\include', ...
            '/usr/local/include' ...
            '/usr/local/mbdyn/include', ...
            '/opt/mbdyn/include', ...
            '/opt/include' ...
                         };

        comparch = computer ('arch');
    
    else
        
        libdir_candidate_locs = {};

        includedir_candidate_locs = {};

        comparch = options.ForceLocalDirsOnlyForArch;
        
    end
    
    switch comparch

        case {'win64', 'mingw32-x86_64'}
            libdir_candidate_locs = [libdir_candidate_locs, ...
                fullfile(getmfilepath ('mexmbdyn_setup'), 'x86_64-w64-mingw32', 'lib')];
            includedir_candidate_locs = [includedir_candidate_locs, ...
                fullfile(getmfilepath ('mexmbdyn_setup'), 'x86_64-w64-mingw32', 'include')];
        case 'win32'
            libdir_candidate_locs = [libdir_candidate_locs, ...
                fullfile(getmfilepath ('mexmbdyn_setup'), 'i686-w64-mingw32', 'lib')];
            includedir_candidate_locs = [includedir_candidate_locs, ...
                fullfile(getmfilepath ('mexmbdyn_setup'), 'i686-w64-mingw32', 'include')];
        case 'glnxa64'
            libdir_candidate_locs = [libdir_candidate_locs, ...
                fullfile(getmfilepath ('mexmbdyn_setup'), 'x86_64-linux-gnu', 'lib')];
            includedir_candidate_locs = [includedir_candidate_locs, ...
                fullfile(getmfilepath ('mexmbdyn_setup'), 'x86_64-linux-gnu', 'include')];
        case 'glnxa32'
            libdir_candidate_locs = [libdir_candidate_locs, ...
                fullfile(getmfilepath ('mexmbdyn_setup'), 'i686-linux-gnu', 'lib')];
            includedir_candidate_locs = [includedir_candidate_locs, ...
                fullfile(getmfilepath ('mexmbdyn_setup'), 'i686-linux-gnu', 'include')];

    end

    if exist ('rnfoundry_setup', 'file') == 2
        switch comparch

            case {'win64', 'mingw32-x86_64'}
                libdir_candidate_locs = [libdir_candidate_locs, ...
                    fullfile(getmfilepath ('rnfoundry_setup'), 'x86_64-w64-mingw32', 'lib')];
                includedir_candidate_locs = [includedir_candidate_locs, ...
                    fullfile(getmfilepath ('rnfoundry_setup'), 'x86_64-w64-mingw32', 'include')];
            case 'win32'
                libdir_candidate_locs = [libdir_candidate_locs, ...
                    fullfile(getmfilepath ('rnfoundry_setup'), 'i686-w64-mingw32', 'lib')];
                includedir_candidate_locs = [includedir_candidate_locs, ...
                    fullfile(getmfilepath ('rnfoundry_setup'), 'i686-w64-mingw32', 'include')];
            case 'glnxa64'
                libdir_candidate_locs = [libdir_candidate_locs, ...
                    fullfile(getmfilepath ('rnfoundry_setup'), 'x86_64-linux-gnu', 'lib')];
                includedir_candidate_locs = [includedir_candidate_locs, ...
                    fullfile(getmfilepath ('rnfoundry_setup'), 'x86_64-linux-gnu', 'include')];
            case 'glnxa32'
                libdir_candidate_locs = [libdir_candidate_locs, ...
                    fullfile(getmfilepath ('rnfoundry_setup'), 'i686-linux-gnu', 'lib')];
                includedir_candidate_locs = [includedir_candidate_locs, ...
                    fullfile(getmfilepath ('rnfoundry_setup'), 'i686-linux-gnu', 'include')];

        end
    end

    libwasfound = false;
    
    libdir = '';
    for ind = 1:numel (libdir_candidate_locs)
        
        if exist (libdir_candidate_locs{ind}, 'dir') == 7
            
            if ( exist (fullfile (libdir_candidate_locs{ind}, 'libmbc.so'), 'file') == 2 ) ...
                || ( exist (fullfile (libdir_candidate_locs{ind}, 'libmbc.a'), 'file') == 2 )

                libdir = libdir_candidate_locs{ind};
                libwasfound = true;
                break;
                
            end
            
        end
        
    end
    
    headerwasfound = false;
    
    includedir = '';
    for ind = 1:numel (includedir_candidate_locs)
        
        if exist (includedir_candidate_locs{ind}, 'dir') == 7
            
            if ( exist (fullfile (includedir_candidate_locs{ind}, 'mbc.h'), 'file') == 2 ) ...
                || ( exist (fullfile (includedir_candidate_locs{ind}, 'mbcxx.h'), 'file') == 2 )

                includedir = includedir_candidate_locs{ind};
                headerwasfound = true;
                break;
                
            end
            
        end
        
    end

end