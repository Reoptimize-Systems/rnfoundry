function rnfoundry_release (varargin)
% perform tasks necessary for a release of the renewnet foundry

    options.W64CrossBuildMexLibsDir = fullfile ('/home', 'rcrozier', 'Sync', 'work', 'matlab_windows_libs', ...
                                               [ 'r', version('-release') ], ...
                                               'extern', 'lib', 'win64', 'mingw64');
	options.RunTests = false;
    options.Version = 'devel';
    options.MBDynVersion = 'devel';
    options.EWSTVersion = 'devel';
    options.ThrowBuildErrors = true;
    options.Verbose = false;
                                           
	options = parse_pv_pairs (options, varargin);
    
    thisfilepath = fileparts (which ('rnfoundry_setup'));
    
    addpath (thisfilepath);
    
    % first build the linux mex fles etc.
    rnfoundry_setup ( 'ForceAllMex', true, ...
                      'Verbose', true, ...
                      'RunTests', options.RunTests, ...
                      'PreventXFemmCheck', true, ...
                      'ThrowBuildErrors', true, ...
                      'Verbose', options.Verbose );
    
	
    % now build the windows mex fles etc.
    [ MBCLibDir, MBCIncludeDir ] = mbdyn.mint.find_libmbc ('ForceLocalDirsOnlyForArch', 'win64');
    
    rnfoundry_setup ( 'ForceAllMex', true, ...
                      'W64CrossBuild', true, ...
                      'W64CrossBuildMexLibsDir', options.W64CrossBuildMexLibsDir, ...
                      'MBCLibDir', MBCLibDir, ...
                      'MBCIncludeDir', MBCIncludeDir, ...
                      'PreventXFemmCheck', true, ...
                      'ThrowBuildErrors', true, ...
                      'Verbose', options.Verbose );
                  
	exampledirs = {};
    
	% now build the documentation
    
    rndocdir = fullfile (thisfilepath, 'documentation');
    
    mkdir (rndocdir);
    
    [ mbdyndocrootdir, mbdyndoczipfilename ] = mbdyn.makedocs ('Version', options.MBDynVersion);
    
    movefile ( fullfile (mbdyndocrootdir, '..', mbdyndoczipfilename), ...
               rndocdir );
           
	[ ewstdocrootdir, ewstdoczipfilename ] = wsim.makedocs ('Version', options.EWSTVersion);
    
    movefile ( fullfile (ewstdocrootdir, '..', ewstdoczipfilename), ...
               rndocdir );
    
    % create the readme files
    help2txtfile (fullfile (thisfilepath, 'README.txt'), 'rnfoundry_release>readme');
    help2txtfile (fullfile (thisfilepath, 'README.rst'), 'rnfoundry_release>readme');

end


function readme ()
% RENEWNET FOUNRDY
% ****************
% 
% The renewnet foundry is a set of matlab codes and other software
% aimed at the simulation of renewable energy systems, including
% development of the power take-off components. The most well
% developed parts of the foundry are the multibody dynamics tools,
% wave systems tools and the permanent magnet machine design and
% simulation tools.
% 
% Most code in the foundry is also compatible with Octave, the free
% alternative to Matlab.
% 
% INSTALLATION
% ============
% 
% To get started with the matlab code you should run the function::
% 
%   rnfoundry_setup
% 
% found in the top level directory. It is also advised that you
% READ THE HELP FOR THIS FUNCTION before running it to get an idea
% of what it does and what the system requirements are to get
% optimum performance. Ideally you will have a C++ compiler set up
% with your matlab installation, see the rnfoundry_setup help for
% more information. Other than this the setup is fully automated.
% 
% USAGE
% =====
% 
% You can find some example scripts in the following directories:
% 
% Edinburgh Wave Systems Toolbox
% ------------------------------
% 
% rnfoundry/wave/doc/sphinx/examples
% 
% This contains examples of using the Edinburgh Wave Systems Toolbox.
% Further (less refined and well commented) examples of using this may
% be found in
% 
% rnfoundry/wave/matlab-octave/wec-sim/test
% 
% 
% MBdyn Multibody Dynamics Toolbox
% --------------------------------
% 
% rnfoundry/common/multibody/MBDyn/doc/sphinx/examples
% 
% Further (less refined and commented) examples can be found in the
% testing code in
% 
% rnfoundry/common/multibody/MBDyn/test
% 
% 
% Permanent Magnet Machines Toolbox
% ---------------------------------
% 
% rnfoundry/common/electrical/matlab-octave/permanent_magnet_machines_tools/examples_and_tutorials
% 
% Contains examples of using the permanent magnet machine simulation
% tools


end


function params = parse_pv_pairs(params,pv_pairs)
% parse_pv_pairs: parses sets of property value pairs, allows defaults
% usage: params=parse_pv_pairs(default_params,pv_pairs)
%
% arguments: (input)
%  default_params - structure, with one field for every potential
%             property/value pair. Each field will contain the default
%             value for that property. If no default is supplied for a
%             given property, then that field must be empty.
%
%  pv_array - cell array of property/value pairs.
%             Case is ignored when comparing properties to the list
%             of field names. Also, any unambiguous shortening of a
%             field/property name is allowed.
%
% arguments: (output)
%  params   - parameter struct that reflects any updated property/value
%             pairs in the pv_array.
%
% Example usage:
% First, set default values for the parameters. Assume we
% have four parameters that we wish to use optionally in
% the function examplefun.
%
%  - 'viscosity', which will have a default value of 1
%  - 'volume', which will default to 1
%  - 'pie' - which will have default value 3.141592653589793
%  - 'description' - a text field, left empty by default
%
% The first argument to examplefun is one which will always be
% supplied.
%
%   function examplefun(dummyarg1,varargin)
%   params.Viscosity = 1;
%   params.Volume = 1;
%   params.Pie = 3.141592653589793
%
%   params.Description = '';
%   params=parse_pv_pairs(params,varargin);
%   params
%
% Use examplefun, overriding the defaults for 'pie', 'viscosity'
% and 'description'. The 'volume' parameter is left at its default.
%
%   examplefun(rand(10),'vis',10,'pie',3,'Description','Hello world')
%
% params = 
%     Viscosity: 10
%        Volume: 1
%           Pie: 3
%   Description: 'Hello world'
%
% Note that capitalization was ignored, and the property 'viscosity'
% was truncated as supplied. Also note that the order the pairs were
% supplied was arbitrary.

    npv = length(pv_pairs);
    n = npv/2;

    if n~=floor(n)
      error 'Property/value pairs must come in PAIRS.'
    end
    if n<=0
      % just return the defaults
      return
    end

    if ~isstruct(params)
      error 'No structure for defaults was supplied'
    end

    % there was at least one pv pair. process any supplied
    propnames = fieldnames(params);
    lpropnames = lower(propnames);
    for i=1:n
      p_i = lower(pv_pairs{2*i-1});
      v_i = pv_pairs{2*i};

      ind = strmatch(p_i,lpropnames,'exact');
      if isempty(ind)
        ind = find(strncmp(p_i,lpropnames,length(p_i)));
        if isempty(ind)
          error(['No matching property found for: ',pv_pairs{2*i-1}])
        elseif length(ind)>1
          error(['Ambiguous property name: ',pv_pairs{2*i-1}])
        end
      end
      p_i = propnames{ind};

      % override the corresponding default in params
      params = setfield(params,p_i,v_i);

    end

end



