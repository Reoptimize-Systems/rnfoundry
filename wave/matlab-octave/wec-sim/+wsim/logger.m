classdef logger < handle
% class to log data with variables of any size
%
% Description
%
% A matlab class for logging data and/or various events during the
% execution of a matlab script/function. This class is mainly a container
% class for logged data, with some additional functionality to operate on
% the stored data in generating plots, applying functions etc.
% 
% The main aim of this class is to consolidate and store all the outputs
% and messages from a MATLAB script/function into a single object.  These
% objects can be saved, and retrieved later, with all the information in
% one place Several utility functions (mainly for plotting) make it easy to
% operate on the stored data.
%
% There are several ways in which you can use this class. You can either
% (1) create a logger object, and start logging into the class
% (2) user's class can be inherited from logger
% (3) you can use matlab's event framework to log events by adding 
% appropriate listeners.
% 
% wsim.logger Methods:
%
%  logger - constructs a logger object
%  addVariable - adds a variable to be logged (initialises data structures)
%  getInfo - get information in a structure for one or more logged variables
%  logErr - log an error  message
%  logMesg - log a general message
%  logWarn - log a warning message
%  logObj - append and object to a logged object series
%  logVal - append a value to a logged variable to the 
%  plot2Vars - plot two scalar logged variables against each other
%  plotFofVar - plot funcitno of variable
%  plotVar - plot a scalar variable
%  setDefaultDesc - set the default description
%  setDesc - set the description for a logged variable
%  setMesgFunc - set the function called when logging a message
%  setPlotFunc - set the function to be used for plotting
%  setSeries - set an entire data series directly
%  setSilent - turn off all messages except errors and warnings
%  truncateVariable - strip preallocated space from variable
%  truncateAllVariables - strip preallocated space from all variables
% 
% Examples
%
%
% Example 1
%
%
% lg = wsim.logger;
% 
% lg.addVariable ( 'weight', [1,1], ...
%                  'Description', 'Weight of Subjects', ...
%                  'AxisLabel', 'kg' );
%              
% lg.addVariable ( 'height', [1,1], ...
%                  'Description', 'Weight of Subjects', ...
%                  'AxisLabel', 'm' );
%              
% 
% for i = 1:100
%     
%     my_output_1 = 10*rand;
%     height = 1.5*my_output_1 + 5*rand;
%     
%     lg.logVal('weight', my_output_1);
%     lg.logVal('height', height);
%     
% end
% 
% 
% lg.plotVar('weight',  'PlotFcnArgs', {'LineWidth', 2, 'Color','r'});
% lg.plotVar('height',  'PlotFcnArgs', {'LineWidth', 2, 'Color','r'});
% lg.plot2Vars('weight', 'height', 'PlotFcnArgs', {'LineWidth', 2, 'Color','r'});
% lg.plotFofVar('height', @log,  'PlotFcnArgs', {'LineWidth', 2, 'Color','r'});
% 
% vals = lg.lastLoggedVals ('weight', 5)
%
%
% Example 2
%
%
% lg = wsim.logger;
% 
% lg.addVariable ( 'weight', [1,1], ...
%                  'Description', 'Weight of Subjects', ...
%                  'AxisLabel', 'kg', ...
%                  'Windowed', true, ...
%                  'PreallocateStorage', 20 );
%              
% lg.addVariable ( 'height', [1,1], ...
%                  'Description', 'Weight of Subjects', ...
%                  'AxisLabel', 'm', ...
%                  'Windowed', true, ...
%                  'PreallocateStorage', 20  );
%              
% 
% for i = 1:100
%     
%     my_output_1 = 10*rand;
%     height = 1.5*my_output_1 + 5*rand;
%     
%     lg.logVal('weight', my_output_1);
%     lg.logVal('height', height);
%     
% end
% 
% lg.plot2Vars('weight','height', 'PlotFcnArgs', {'LineWidth', 2, 'Color','r'});
% lg.plotVar('weight',  'PlotFcnArgs', {'LineWidth', 2, 'Color','r'});
% lg.plotVar('height',  'PlotFcnArgs', {'LineWidth', 2, 'Color','r'});
% lg.plotFofVar('height',@log,  'PlotFcnArgs', {'LineWidth', 2, 'Color','r'});
% 
% % this will return only 20 values, because that's all there is
% vals = lg.lastLoggedVals ('weight', 25)
% 
%
% Author: Richard Crozier
%
% Derived from the original logger class developed by Pavan Mallapragada 
% Organization:  Massachusetts Institute of Technology
% Contact:       <pavan_m@mit.edu>
%


    properties 
        data; % A structure in which the logged data will be stored.
    end
    
    properties (GetAccess = public, SetAccess = private)
        
        % A structure that stores the descriptions of the fields. These descriptions are used for labeling plots.
        info;
        
        % Number of variables being logged
        numVariables;

        % A cell array that stores all the fields that are being logged by the object currently, as strings.
        fieldNames;

        % A string containing field name to perform default operations on.  Currently this is unused.
        defaultfield;

        % Cell array to store messages
        messages;

        % Cell array to store warnings
        warnings;

        % Cell array to store error messages
        errors;

        % Boolean value. If set to true, logger will not print any of its messages. It will still print the warnings/errors that the user specifies it to.
        silent;

        % The function handle used for plotting. Default is plot.
        plotfunc;

        % The function handle used for displaying interal messages of logger object. Default is set to warning.
        mesgfunc;

        % A string used for labeling the x-axis when plotting single variables.
        defaultXLabel;
        
        % factor by which to expand preallocated space when the end of the
        % preallocated space is reached and more data is logged
        expandPreallocFactor;
        
    end
    

    methods 

        function obj = logger(varargin)
            % constructs a logger object.
            %
            % Syntax
            %
            % obj = wsim.logger ()
            % obj = wsim.logger ('Parameter', Value)
            %
            % Description
            %
            % Initializes basic plotting functions, verbosity and messaging
            % options. Initializes the structures used for storing
            % messages/warnings/errors etc.
            %
            % Input
            %
            % Arguments may be supplied as parameter-value pairs. The
            % available options are:
            %
            %  'ExpandPreallocFactor' - variables can have data
            %    preallocated to improve speed. This option can be used to
            %    control how the preallocated arrays expand when the end of
            %    the preallocated array is reached. When this occurs the
            %    arrays are expanded by this factor multiplied by the
            %    current array length. Default is 0.25 if not supplied.
            %
            %  'PlotFcn' - handle to function to be used for plotting
            %    variables. Default is @plot if not supplied.
            %
            %  'MesgFcn' - handle to function to be called when reporting
            %    messages. By default this is @warning, but it could, for
            %    example, be made  @error to result in more stringent
            %    checking.
            %
            %  'Silent' - Sets the silent property true or false. If true,
            %    then the internal error messages of the logger class are
            %    not printed out. Whatever the user specifies by setting
            %    the show variable in the logWarn, logMesg functions are
            %    still printed out to the command window. Defualt is true.
            %
            %  'DefaultDescription' - character vector containing the
            %    default description for variables, which is used as the x
            %    axis label in plots. Default is '' (an empty character
            %    vector).
            %
            % Output
            %
            %  obj - wsim.logger object
            %
            %

            
            options.ExpandPreallocFactor = 0.25;
            options.PlotFcn = @plot;
            options.MesgFcn = @warning;
            options.Silent = true;
            options.DefaultDescription = '';
            
            options = parse_pv_pairs (options, varargin);
            
            check.isNumericScalar (options.ExpandPreallocFactor, true, 'ExpandPreallocFactor', 1);
            assert (isa (options.PlotFcn, 'function_handle'), ...
                'PlotFcn must be a function handle' );
            assert (isa (options.MesgFcn, 'function_handle'), ...
                'MesgFcn must be a function handle' );
            check.isLogicalScalar (options.Silent, true, 'Silent');
            assert (ischar (options.DefaultDescription), ...
                'DefaultDescription must be a character vector' );
            
            obj.silent = false;

            obj.info = struct ();

            obj.messages = {};
            obj.warnings = {};
            obj.errors   = {};

            obj.plotfunc = options.PlotFcn;
            obj.mesgfunc = options.MesgFcn;
            obj.silent   = options.Silent;

            obj.defaultXLabel = options.DefaultDescription;
            
            obj.numVariables = 0;
            obj.fieldNames = {};
            
        end
        
        
        function addVariable (obj, name, varsize, varargin)
            % add a variable to be logged
            %
            % Syntax
            %
            % addVariable (lo, name, varsize)
            % addVariable (..., 'Parameter', value)
            %
            % Description
            %
            % addVariable adds a variable to the logger. Variables can
            % scalars, vectors or matrices of any size. In the case of
            % matrices and vectors logging is performed by expanding the
            % variables along the first singleton dimension (which will be
            % created if there are no singleton dimensions). Storage for
            % the logged data can be preallocated to improve speed.
            %
            % Input
            %
            %  lo - wsim.logger object
            %
            %  name - name of the variable to be logged
            %
            %  varsize - size of the variable to be logged in the same
            %   format as would be reported by the 'size' function, i.e. a
            %   vector of integers. The log of the variable will be made by
            %   expanding the first singleton dimension (dimension of size
            %   1) of the variable, i.e. new values of the variable at each
            %   logging step will be appended to this dimension. If there
            %   are no singleton dimensions, an extra dimension will be
            %   created. The logging dimesion can also be specified using
            %   the 'ForceLogDimension' option described below.
            %
            % Additional arguments can be supplied using parameter-value
            % pairs. The available options are:
            %
            %  'Description' - optional string containing a description of
            %    the variable
            %
            %  'ForceLogDimension' - optional scalar integer giving the
            %    dimension along which the logged variable is to be
            %    appended at each logging step.
            %
            %  'PreallocateStorage' - optional scalar integer giving the
            %    number of samples for which to preallocate storage in
            %    advance. When this limit is reached, the storage space is
            %    expanded, unless the 'Windowed' option is used (see
            %    below).
            %
            %  'IndependentVariable' - string containing the name of
            %    another variable, already added to the logger, which is
            %    the independant variable for the new variable. Primarily
            %    used for plotting, but, also for preallocation. If the
            %    PreallocateStorage option is not specified, the new
            %    variable is by default preallocated the same storage as
            %    the independant variable.
            %
            %  'AxisLabel' - optional charachter vector which will be used
            %    as the axis label for this variable when plotting it. Will
            %    appear on the axis against which this variable is plotted,
            %    which could be either X or Y.
            %
            %  'Windowed' - optional true/false flag indicating whether the
            %    data should be 'windowed' i.e. once the number of values
            %    logged equals the value in PreallocateStorage (see above),
            %    instead of increasing storage space, the earlier values
            %    will be overwritten with new data, where the earliest
            %    logged value is lost with each new value, and all current
            %    values are shifted earlier in the log. Default is false,
            %    meaning the storage space will increase once the
            %    preallocated space is reached.
            %
            %

            
            options.Description = '';
            options.AxisLabel = '';
            options.ForceLogDimension = [];
            options.PreallocateStorage = [];
            options.IndependentVariable = '';
            options.Windowed = false;
            options.Legends = {};
            
            options = parse_pv_pairs (options, varargin);
            
            check.isLogicalScalar (options.Windowed, true, 'Windowed');
            
            assert (isvarname (name), 'name must be a string containing a valid variable name');
            
            assert (isvector (varsize) && numel (varsize) >= 2, ...
                'varsize must be a vector of variable dimensions with two or more elements' );
            
            if ~isempty (options.PreallocateStorage)
                check.isPositiveScalarInteger (options.PreallocateStorage, true, 'PreallocateStorage');
            end
            
            assert (ischar (options.Description), 'Description must be a character vector');
            assert (ischar (options.AxisLabel), 'AxisLabel must be a character vector');

            if ~isempty (options.IndependentVariable)
                
                assert (isvarname (options.IndependentVariable), ...
                    'IndependentVariable must be a string containing a valid variable name');
                
                assert (isfield (obj.info, options.IndependentVariable), [...
'If specifiying an independent variable for a logged variable, you must add\n' ...
'the independent variable first. The specified independent variable, "%s", \n' ...
'was not found in this logger object.' ]);

                if isempty (options.PreallocateStorage)
                    % by default preallocate the same storeage as the
                    % independant variable if the preallocation is not
                    % explicitly specified
                    options.PreallocateStorage = obj.info.(options.IndependentVariable).PreallocatedLogLength;
                    
                end
                
            else
                
                if isempty (options.PreallocateStorage)
                    % default preallocation is 1
                    options.PreallocateStorage = 1;
                end
                
            end
            
            % generate input structure for subasgn function, adding
            % to log is done using this function to add data along the
            % first dimesion with length 1
            indass.type = '()';
            logdim = 0;
            loggedvarsize = varsize;
            datadims = 1:numel(varsize);
            
            if ~isempty (options.ForceLogDimension)
                
                check.isPositiveScalarInteger (options.ForceLogDimension, true, 'ForceLogDimension');
                
                % use user specified dimension
                if options.ForceLogDimension > numel (loggedvarsize)
                    % add the necessary singleton dimensions to varsize
                    for ind = (numel (loggedvarsize) + 1):options.ForceLogDimension
                        loggedvarsize(ind) = 1;
                    end
                else
                    % check this dimension is of size 1
                    assert (loggedvarsize(options.ForceLogDimension) == 1, ...
                        'The size of data dimension chosen in ForceLogDimension must be 1');
                end
                
                indass.subs = repmat ({':'}, 1, numel (loggedvarsize));
                indass.subs{options.ForceLogDimension} = 1;
                logdim = options.ForceLogDimension;
                    
            else
                % find the first singleton dimension of array and log along
                % it
                
                for ind = 1:numel (loggedvarsize)
                    if loggedvarsize(ind) == 1
                        logdim = ind;
                    end
                end
                
                if logdim == 0
                    % no singleton dimensions were found, so add one at the
                    % end
                    loggedvarsize(end+1) = 1;
                    logdim = numel (loggedvarsize);
                    indass.subs = repmat ({':'}, 1, numel (loggedvarsize));
                else
                    if numel (loggedvarsize) == 2 && loggedvarsize(1) == 1
                        indass.subs = {1, ':'};
                        logdim = 1;
                        datadims = 2;
                    elseif numel (loggedvarsize) == 2 && loggedvarsize(2) == 1
                        indass.subs = {':', 1};
                        logdim = 2;
                        datadims = 1;
                    else
                        indass.subs = repmat ({':'}, 1, numel (loggedvarsize));
                    end
                end
                
                indass.subs{logdim} = 1;
            
            end
            
            if isfield (obj.info, name)
                error ('A variable named %s is already being logged, you must choose a unique name.', name);
            else
                
                obj.info.(name) = struct ( 'Description', options.Description, ...
                                           'Size', varsize, ...
                                           'PreallocatedLogLength', options.PreallocateStorage, ...
                                           'IndexAssignment', indass, ...
                                           'IndexDimension', logdim, ...
                                           'DataDimensions', datadims, ...
                                           'LastLogIndex', 0, ...
                                           'IndependentVariable', options.IndependentVariable, ...
                                           'LoggedVariableNumber', numel(obj.info) + 1, ...
                                           'LoggedSize', loggedvarsize, ...
                                           'Windowed', options.Windowed, ...
                                           'AxisLabel', options.AxisLabel, ...
                                           'Legends', {options.Legends} );

            end
            
            % make the singleton logging dimension the required
            % length, the data will be preallocated with nans
            loggedvarsize(logdim) = options.PreallocateStorage;

            obj.data.(name) = nan (loggedvarsize);
            
            obj.fieldNames = fieldnames (obj.info);
            obj.numVariables = numel (obj.fieldNames);
            
        end
 

        function [] = logErr(obj, mesg, show)
            % log an error text message 
            %
            % Syntax
            %
            % logErr(obj, mesg, show)
            %
            % Description
            %
            % Logs text messages given as an error, which are stored in
            % the errors property of the logger class.
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  mesg - character vector to be stored in the errors log.
            %
            %  show - true/false flag indicating whether to also print the
            %   message to the command line.
            %
            %
            % See Also: wsim.logger.logWarn, wsim.logger.logMesg
            %
            
            assert (ischar (mesg), 'mesg must be a character vector');
            
            obj.errors{end+1} = mesg;
            if nargin > 2 && show
                fprintf('(e): %s\n',mesg);
            end
        end
        
 
        function [] = logWarn(obj, mesg, show)
            % log a warning text message 
            %
            % Syntax
            %
            % logWarn(obj, mesg, show)
            %
            % Description
            %
            % Logs text messages given as a warning, which are stored in
            % the warnings property of the logger class.
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  mesg - character vector to be stored in the warnings log.
            %
            %  show - true/false flag indicating whether to also print the
            %   message to the command line.
            %
            %
            % See Also: wsim.logger.logErr, wsim.logger.logMesg
            %
            
            assert (ischar (mesg), 'mesg must be a character vector');
            
            obj.warnings{end+1} = mesg;
            if nargin > 2 && show
                fprintf('(w): %s\n',mesg);
            end
        end
        
 
        function logMesg(obj, mesg, show)
            % log a general text message 
            %
            % Syntax
            %
            % logMesg(obj, mesg, show)
            %
            % Description
            %
            % Logs text messages given for information, which are stored in
            % the messages property of the logger class.
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  mesg - character vector to be stored in the message log.
            %
            %  show - true/false flag indicating whether to also print the
            %   message to the command line.
            %
            %
            % See Also: wsim.logger.logErr, wsim.logger.logWarn
            %
            
            assert (ischar (mesg), 'mesg must be a character vector');

            obj.messages{end+1} = mesg;
            if nargin > 2 && show
                fprintf('(i): %s\n', mesg);
            end
            
        end

        
        function [] = logObj(obj, field, obj2log)
            % logs non-numeric objects, and stores them in a field as a cell array.
            %
            % Description
            %
            % Log a non-numeric object. Non-numeric object data cannot be
            % plotted.
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  field - The field name you are logging.
            %
            %  val - The value any matlab object.
    
            try
                obj.data.(field){end+1} = obj2log;
            catch
                obj.addprop(field);
                if ~obj.silent
                    fprintf(['Logger: new field added to the logger object: ', field,'\n']);
                end
                obj.fieldNames{end+1} = field;
                obj.data.(field) = {};
                obj.data.(field){end+1} = obj2log;
            end

        end


        function setMesgFunc(obj, mf) 
            % set the default waring/error function for the logger object
            %
            % Syntax
            %
            % setMesgFunc(obj, mf)
            %
            % Description
            %
            % setMesgFunc sets the function which is run when an internal
            % warning is issued (this is NOT related to the logMesg,
            % logWarn or logErr methods). By default @warning is used for
            % most internal warnings/errors, but this can be replaced, e.g.
            % with @error for more strict treatment of problems.
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  mf - function handle to new function to be run when an
            %   internal error occurs. Should take a single character
            %   vector input.
            %
            %
            
            assert (isa (mf, 'function_handle'), ...
                'mf must be a function handle');
            
            obj.mesgfunc  = mf;
            
        end 


        function [] = setDefaultDesc(obj, str)
            % set the default description which is used for labeling the x-axis.
            %
            % Description
            %
            % The default description is used when plotting, for the X axis
            % label when there is no independant variable associated with a
            % variable.
            %
            
            
            obj.defaultXLabel = str;
        end


        function [] = disp(obj)
            % Overriden disp function for the logger class.
            %
            % Description
            %
            % Displays the logger object with all its variables and their
            % sizes.
            %
            
            fprintf('Logger Object\n\n');
            fprintf('Logged User Specified Variables:\n');

            wid = max (cellfun (@length, obj.fieldNames, 'UniformOutput', true)) + 4;
            
            printstr = sprintf ('%%%ds : %%d log entries\\n', wid);
            
            
            for i = 1:numel(obj.fieldNames)
                fprintf(1, ...
                        printstr, ...
                        obj.fieldNames{i}, ...
                        obj.info.(obj.fieldNames{i}).LastLogIndex );
            end

            fprintf('\nOther Logged Data:\n');
            

            fprintf(1, printstr, 'messages', length(obj.messages));
            fprintf(1, printstr, 'warnings', length(obj.warnings));
            fprintf(1, printstr, 'errors', length(obj.errors));
        end


        function setDesc(obj, varname, desc)
            % set the description for a variable to be used for x-y labels in plots
            %
            % Syntax
            %
            % setDesc(obj, varname, desc)
            %
            % Description
            %
            % setDesc set the description for one or more variables. The
            % description is used when plotting as the label for the
            % appropriate axis.
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  varname - either a character vector containing the name of
            %   logged variable for which the description is to be set, or
            %   a cell array of character vectors containing the names of
            %   several variables for which the descriptions are to be set.
            %   If a cell array, it must be the same size as desc (see
            %   below).
            %
            %  desc - character vector containing the new description for
            %   the variable in varname, or a cell array of character
            %   vectors containing a description for each variable in
            %   varname (when it is a cell array).
            %

            if ischar(varname) && ischar(desc)
                obj.info.(varname).Description = desc;
            elseif iscell(varname) && iscell(desc)
                
                assert (samesize (varname, desc), ...
                    'If varname and desc are cell string arrays they must be the same size');
                
                for i = 1:numel(varname)
                    obj.info.(varname{i}).Description = desc{i};
                end
            else
                error ('varname and desc must be both either character vectors, or cell arrays of the same size')
            end
            
        end


        function setPlotFunc(obj, pf_handle)
            % Set the plotting function to be used. 
            %
            %
            % Syntax
            %
            % setPlotFunc(obj, pf_handle)
            %
            % Description
            %
            % Set the plotting function to be used when calling plotVar,
            % plotFofVar, plot2Vars etc.
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  pf_handle - function handle containing plotting function to
            %   be used.
            %
            % See Also: wsim.logger.plotVar, wsim.logger.plotFofVar
            %           wsimlogger.plot2Vars
            %
            
            obj.plotfunc = pf_handle;
        end


        function [hplot, hax, hfig] = plotFofVar (obj, varname, func, varargin)
            % apply function to logged variable and plot the result
            %
            %
            % Syntax
            %
            % h = plotFofVar(obj, varname, func) 
            % h = plotFofVar(obj, varname, func, plotarg1, plotarg2, ..., plotargn)
            %
            % Description
            %
            % plotFofVar applies a function to a logged variable and plots
            % the result. For example, if a user logs a variable 'height',
            % the log(heights) can be plotted as
            % logger.plotFofVar('height',@log). The parameters to the plot
            % can be provided after the fields, and all those arguments go
            % to the plot function.  For example,
            % logger.plotFofVar('height','@log','LineWidth',2,'Color','r');
            % will pass the last four arguments to plot.
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  varname - name of logged variable to be plotted
            %
            %  func - function handle containing function to be applied to
            %   the variable before plotting.
            %
            % Output
            %
            %  hplot - A handle (or array of handles) for the plot series
            %   generated. Useful for formatting by the user.
            %
            %  hax - handle to the axes object in which the plot was drawn.
            %
            %  hfig - handle to the figure object in which the plot was
            %   drawn
            %
            %
            % See Also: wsim.logger.plotVar, wsim.logger.plot2Vars
            %
            
            options.PlotFcnArgs = {};
            options.Skip = 1;
            
            options = parse_pv_pairs (options, varargin);

            if ~ischar(varname), obj.mesgfunc([varname 'must be a string specifying field that are already added to the logger object.']); return; end
            if ~isnumeric(obj.data.(varname)) obj.mesgfunc(['Plotting only numeric values is supported at this point. Not generating the plot for' varname]); return; end
            if ~isa(func, 'function_handle') obj.mesgfunc('third argument func must be a function handle'); return; end

            plotvals = feval (func, obj.data.(varname));
            
            hfig = figure ();
            hax = axes ();

            hplot = obj.plotfunc (hax, obj.data.(varname), options.PlotFcnArgs{:});

            hold on;

            ylabel ([func2str(func) ' (' obj.info.(varname).AxisLabel ')'], 'FontSize', 16);
            if ~isempty (obj.defaultXLabel)
                xlabel (obj.defaultXLabel,'FontSize',16);
            end

            set (hax, 'FontSize', 16);

            hold off;
            
        end


        function [hplot, hax, hfig] = plotVar (obj, varname, varargin)
            % given a string specifying a numeric scalar field, it is plotted.
            %
            % Syntax
            %
            % plotVar(obj, field)
            % plotVar(obj, field, plotarg1, plotarg2, ..., plotargn)
            % [hplot, hax, hfig] = plotVar(...)
            %
            % Description
            %
            % plotVar Plots one of the logged variables. If a user logs a
            % fields 'height', using logger.plotVar('height') will plot
            % height. Additional arguments for the plotting funtion (e.g.
            % plot) can be provided using the 'PlotFcnArgs' option.
            %
            % Input
            %
            %  obj - logger object
            % 
            %  field - a string specifying the name of logged field 1.
            %
            % Additional arguments may be supplied as parameter-value
            % pairs. The available options are:
            %
            %  'PlotFcnArgs' - cell array of additional arguments that are 
            %    to be sent to the plotting function.
            %
            %  'Skip' - integer defining a number of values to skip between
            %    plotted values (to reduce plot size for large series)
            %
            %  'Scale' - scalar value by which to scale the data values
            %    when plotting
            %
            % Output
            %
            %  hplot - A handle (or array of handles) for the plot series
            %   generated. Useful for formatting by the user.
            %
            %  hax - handle to the axes object in which the plot was drawn.
            %
            %  hfig - handle to the figure object in which the plot was
            %   drawn
            %
            %
    
            options.PlotFcnArgs = {};
            options.Skip = 1;
            options.Scale = 1;
            
            options = parse_pv_pairs (options, varargin);
            
            if ~ischar(varname), obj.mesgfunc([varname 'must be a string specifying field that are already added to the logger object.']); return; end
            
            if ~isfield (obj.data, varname)
                obj.mesgfunc(['Variable %s does not appear to exist.' varname]);
            end
            
            if ~isnumeric(obj.data.(varname)) obj.mesgfunc(['Plotting only numeric values is supported at this point. Not generating the plot for' varname]); return; end

            indepvar = obj.info.(varname).IndependentVariable;
            
            if isempty (indepvar)
                
                hfig = figure;
                if isoctave
                    hax = axes ();
                else
                    hax = axes (hfig);
                end
                
                hplot = obj.plotfunc (hax, obj.data.(varname) .* options.Scale, options.PlotFcnArgs{:});
            
                if options.Scale == 1
                    ylabel (obj.info.(varname).AxisLabel, 'FontSize', 16);
                else
                    ylabel (sprintf ('%s X %g', obj.info.(varname).AxisLabel, options.Scale), 'FontSize', 16);
                end

                if ~isempty (obj.defaultXLabel)
                    xlabel (obj.defaultXLabel, 'FontSize',16);
                end

                set (hax, 'FontSize', 16);
                
%                 indexvarname = sprintf ('%s_index', varname);
%                 
%                 obj.data = setfield (obj.data, indexvarname, 1:obj.data.(varname).last );
%                 
%                 CC = onCleanup ( @() rmfield (obj.data, indexvarname) );
%                 
%                 
%                 [hplot, hax, hfig] = plot2Vars ( obj, ...
%                                               obj.info.(varname).IndependentVariable, ...
%                                               varname, ...
%                                               'PlotFcnArgs', options.PlotFcnArgs, ...
%                                               'Skip', options.Skip );
%                 
            else
                
                 [hplot, hax, hfig] = plot2Vars ( obj, ...
                                              obj.info.(varname).IndependentVariable, ...
                                              varname, ...
                                              'PlotFcnArgs', options.PlotFcnArgs, ...
                                              'Skip', options.Skip, ...
                                              'Scale', [1, options.Scale]);
                
            end

        end


        function [h, hax, hfig] = plot2Vars (obj, f1, f2, varargin)
            % plot two logged variables against each other
            %
            %
            % Syntax
            %
            % plot2Vars (obj, f1, f2)
            % plot2Vars (obj, f1, f2, plotarg1, plotarg2, ..., plotargn)
            % [h, hax, hfig] = plot2Vars (...)
            %
            % Description
            %
            % plot two logged variables against each other, e.g. if a user
            % logs two variables 'height' and 'weight', giving
            % logger.plotvars('height','weight') will plot height vs.
            % weight. The parameters to the plot can be provided after the
            % variables, and all those arguments go to the plot function.
            % For example,
            %
            % logger.plotvars('height','weight','LineWidth',2,'Color','r');
            %
            % will pass the last four arguments to plot.
            %
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  f1 - name of first variable to plot
            %
            %  f2 - name of second variable to plot
            %
            % Additional arguments may be supplied as parameter-value
            % pairs. The available options are:
            %
            %  'PlotFcnArgs' - cell array of additional arguments that are 
            %    to be sent to the plotting function.
            %
            %  'Skip' - integer defining a number of values to skip between
            %    plotted values (to reduce plot size for large series)
            %
            %  'Scale' - two element vector containing the values by which
            %    to scale the data values for each variable when plotting.
            %    Default is [1,1] if not supplied.
            %
            % Output
            %
            %  h - A handle (or array of handles) for the plot series
            %   generated. Useful for formatting by the user.
            %
            %  hax - handle to the axes object in which the plot was drawn.
            %
            %  hfig - handle to the figure object in which the plot was
            %   drawn
            %
            %
            % See Also: wsim.logger.plotVar
            %
            
            options.PlotFcnArgs = {};
            options.Skip = 1;
            options.Scale = [1,1];
            
            options = parse_pv_pairs (options, varargin);
            
            check.isPositiveScalarInteger (options.Skip, true, 'Skip');

            if ~ischar(f1) || ~ischar(f2)
                feval (obj.mesgfunc, 'f1 and f2 must be strings specifying fields that are already added to the logger object.');
            end
            
            if ~isnumeric(obj.data.(f1)) || ~isnumeric(obj.data.(f2))
                feval (obj.mesgfunc, sprintf ('Plotting only numeric values is supported at this point. Not generating the plot %s vs %s', f1, f2));
            end
            
            f1sz = size (obj.data.(f1));
            f2sz = size (obj.data.(f2));
            
            f1sztest = f1sz;
            f1sztest(obj.info.(f1).IndexDimension) = 1;
            
            if ( numel (obj.info.(f1).Size) ~= numel (obj.info.(f2).Size) ) ...
                    && ( all(f1sztest ~= 1) )
                feval (obj.mesgfunc, 'Cannot plot variables of different sizes against each other (unless one variable is a vector)');
            end
            
            if ndims (obj.data.(f1)) > 3
                feval (obj.mesgfunc, 'Cannot plot variables of more than 2 dimensions');
            end
            
            f2datadims = obj.info.(f2).DataDimensions;
                
            f2szdim1 = size (obj.data.(f2), f2datadims(1));
            
            if numel (f2datadims) > 1
                f2szdim2 = size (obj.data.(f2), f2datadims(2));
            else
                f2szdim2 = 1;
            end
            
            h = [];
            f2S.type = '()';
            f2S.subs = cell (1, numel(f2sz));
            for ind = 1:numel(f2S.subs), f2S.subs{ind} = 1; end
            if options.Skip ~= 1
                f2S.subs{obj.info.(f2).IndexDimension} = 1:options.Skip:obj.info.(f2).LastLogIndex;
            else
                f2S.subs{obj.info.(f2).IndexDimension} = ':';
            end
            
            
            f1datadims = obj.info.(f1).DataDimensions;

            f1S.type = '()';
            f1sz = size (obj.data.(f1));
            f1S.subs = cell (1, numel(f1sz));
            for ind = 1:numel(f1S.subs), f1S.subs{ind} = 1; end
            if options.Skip ~= 1
                f1S.subs{obj.info.(f1).IndexDimension} = 1:options.Skip:obj.info.(f1).LastLogIndex;
            else
                f1S.subs{obj.info.(f1).IndexDimension} = ':';
            end
            
            indevarsamesize = all (obj.info.(f1).Size == obj.info.(f2).Size);
            
            legstrs = {};
            
            % make a new figure to plot in
            hfig = figure;
            if isoctave
                hax = axes ();
            else
                hax = axes (hfig);
            end
            
            hold on
            
            switch numel (f2datadims)
                
                case 1
                    
                    for dataind = 1:size (obj.data.(f2), f2datadims)
                        
                        % get the indexing for this series
                        f2S.subs{f2datadims} = dataind;
                    
                        if indevarsamesize
                            % if the independant variable, f1 is the same size
                            % as the dependant variable, f2, plot each
                            % component of f1 against f2. Otherwise we always
                            % plot the same series from f1 against each
                            % component of f2
                            f1S.subs{f1datadims} = dataind;
                        end

                        h = [ h, obj.plotfunc( hax, ...
                                               squeeze ( subsref (obj.data.(f1), f1S) .* options.Scale(1) ), ...
                                               squeeze ( subsref (obj.data.(f2), f2S) .* options.Scale(2) ), ...
                                               options.PlotFcnArgs{:} ) ...
                            ];
                        
                        legstrs = [ legstrs, {sprintf('Series %d', dataind)}];
                        
                    end
                    
                case 2
                    
                    for dataind1 = 1:size (obj.data.(f2), f2datadims(1))
                        
                        for dataind2 = 1:size (obj.data.(f2), f2datadims(2))
                        
                            % get the indexing for this series
                            f2S.subs{f2datadims(1)} = dataind1;
                            f2S.subs{f2datadims(2)} = dataind2;

                            if indevarsamesize
                                % if the independant variable, f1 is the same size
                                % as the dependant variable, f2, plot each
                                % component of f1 against f2. Otherwise we always
                                % plot the same series from f1 against each
                                % component of f2
                                f1S.subs{f1datadims(1)} = dataind1;
                                f1S.subs{f1datadims(2)} = dataind2;
                            end

                            h = [ h, obj.plotfunc( hax, ...
                                                   squeeze ( subsref (obj.data.(f1), f1S) .* options.Scale(1) ), ...
                                                   squeeze ( subsref (obj.data.(f2), f2S) .* options.Scale(2) ), ...
                                                   options.PlotFcnArgs{:} ) ...
                                ];
                            
                            if isempty (obj.info.(f2).Legends)
                                legstrs = [ legstrs, {sprintf('Series (%d,%d)', dataind1, dataind2)}];
                            else
                                legstrs = [ legstrs, obj.info.(f2).Legends(dataind1, dataind2)];
                            end
                        
                        end
                        
                    end
                    
                    
                otherwise
                    error ('plotting data with more than 2 data dimensions (i.e. three dimensions in total) is not currently supported');
            end
                
%             for datadimind = 1:numel(f2datadims)
%                 
%                 % reset all the non index dimension subs to 1
%                 indexdimsize = f1S.subs{obj.info.(f1).IndexDimension};
%                 for ind = 1:numel(f1S.subs), f1S.subs{ind} = 1; end
%                 f1S.subs{obj.info.(f1).IndexDimension} = indexdimsize;
%                 
%                 indexdimsize = f2S.subs{obj.info.(f2).IndexDimension};
%                 for ind = 1:numel(f2S.subs), f2S.subs{ind} = 1; end
%                 f2S.subs{obj.info.(f2).IndexDimension} = indexdimsize;
%                 
%                 % check if this dimension is the index dimension of the
%                 % dependant variable
%                 if f2datadims(datadimind) ~= obj.info.(f2).IndexDimension
%                     % if it it not the index dimension we have to plot each
%                     % series from this dimension along the index dimension
%                     
%                     for dataind = 1:size (obj.data.(f2), f2datadims(datadimind))
%                         
%                         % get the indexing for this series
%                         f2S.subs{f2datadims(datadimind)} = dataind;
%                     
%                         if indevarsamesize
%                             % if the independant variable, f1 is the same size
%                             % as the dependant variable, f2, plot each
%                             % component of f1 against f2. Otherwise we always
%                             % plot the same series from f1 against each
%                             % component of f2
%                             f1S.subs{f1datadims(datadimind)} = dataind;
%                         end
% 
%                         h = [ h, obj.plotfunc( squeeze ( subsref (obj.data.(f1), f1S) ), ...
%                                                squeeze ( subsref (obj.data.(f2), f2S) ), ...
%                                                options.PlotFcnArgs{:} ) ...
%                             ];
%                         
%                     end
%                     
%                     
%                 end
%                 
%             end

            hold off;
            
            if options.Scale(1) == 1
                desc1 = obj.info.(f1).AxisLabel;
            else
                desc1 = sprintf ('%s X %g', obj.info.(f1).AxisLabel, options.Scale(1));
            end
            
            if options.Scale(2) == 1
                desc2 = obj.info.(f2).AxisLabel;
            else
                desc2 = sprintf ('%s X %g', obj.info.(f2).AxisLabel, options.Scale(2));
            end

            xlabel (desc1, 'FontSize', 12);
            ylabel (desc2, 'FontSize', 12);
            set (hax, 'FontSize', 12);
            
            if numel (legstrs) > 1
                legend (legstrs{:});
            end
            
            
            
        end


        function status = logVal(obj, varname, val, ignoremissing, checkexists)
            % log a new value of a variable
            %
            % Syntax
            %
            % status = logVal(obj, varname, val, ignoremissing)
            %
            % Description
            %
            % If a field exists, the element is added to the next index of
            % the field.
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  varname - name of the variable to be appended
            %
            %  val - value to be appended to the history of values
            %
            %  ignoremissing - optional true/false flag indicating whether
            %   an error should be thrown if varname does not exist in the
            %   logger. If false, an error is thrown if varname does not
            %   exist. If true an error will not be thrown, and the success
            %   is indicated in the status output. Default is false if not
            %   supplied.
            %
            % Output
            %
            %  status - 0 if the variable varname exists, -1 if not.
            %
            % See Also: wsim.logger.addVariable
            %

            status = 0;
            
            if nargin < 4
                ignoremissing = false;
            end
            
            if nargin < 5
                checkexists = true;
            end
            
            if checkexists
                if ~isfield (obj.data, varname)
                    if ignoremissing
                        status = -1;
                        return;
                    else
                        error ('data logging field: %s does not exist', varname);
                    end
                end
            end
            
            S = [];
            
            if obj.info.(varname).LastLogIndex + 1 > obj.info.(varname).PreallocatedLogLength
                
                if obj.info.(varname).Windowed
                    % rotate the log as we have reached the end of the
                    % window
                    obj.data.(varname) = circshift (obj.data.(varname), 1, obj.info.(varname).IndexDimension);
                    % step the LastLogIndex back one, so the new data
                    % replaces the last item in the log (the data is later
                    % placed in the  LastLogIndex + 1 position)
                    obj.info.(varname).LastLogIndex = obj.info.(varname).LastLogIndex - 1;
                else
                    % expand the data array as we have run out of preallocated
                    % space. Use the expandPreallocFactor to determine by how
                    % much to expand
                    obj.info.(varname).PreallocatedLogLength = max ( ...
                                    obj.info.(varname).LastLogIndex + 1, ...
                                    round (obj.info.(varname).LastLogIndex * obj.expandPreallocFactor) ...
                                                                   );
                                                               
                    % copy the pre-constructed indexing structure (created when
                    % adding the variable)
                    S = obj.info.(varname).IndexAssignment;
                    
                    % build the correct index into the logged variable by
                    % replacing the appropriate index with the new index of the
                    % end of the preallocated data
                    S.subs{obj.info.(varname).IndexDimension} = obj.info.(varname).PreallocatedLogLength;

                    % assign nan to expand the array
                    obj.data.(varname) = subsasgn (obj.data.(varname), S, nan);
                end
            
            end
            
           
            
            % assign the new value
            switch obj.info.(varname).IndexDimension
                
                case 1
                    
                    obj.data.(varname)(obj.info.(varname).LastLogIndex + 1) = val;
                    
                case 2
                    
                    obj.data.(varname)(:,obj.info.(varname).LastLogIndex + 1) = val;
                    
                case 3
                    
                    obj.data.(varname)(:,:,obj.info.(varname).LastLogIndex + 1) = val;
                    
                otherwise
                    if isempty (S)
                        % copy the pre-constructed indexing structure (created when
                        % adding the variable)
                        S = obj.info.(varname).IndexAssignment;
                    end
            
                    % build the correct index into the logged variable by replacing
                    % the appropriate index with the new log index
                    S.subs{obj.info.(varname).IndexDimension} = obj.info.(varname).LastLogIndex + 1;
                    
                    obj.data.(varname) = subsasgn (obj.data.(varname), S, val);
                    
            end
            
            % increment the data index counter for this field
            obj.info.(varname).LastLogIndex = obj.info.(varname).LastLogIndex + 1;

        end
        
        
        function vals = lastLoggedVals(obj, varname, n)
            % get last 'n' logged values of a variable
            %
            % Syntax
            %
            % val = lastLoggedVal(obj, varname, n)
            %
            % Description
            %
            % If a variable exists, the last 'n' values logged to that
            % variable are returned. If the number of logged variables is
            % less than this, all available data will be returned.
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  varname - name of the variable to be for which to obtain the
            %    last logged value
            %
            %  n - (optional) the number of logged variables to return.
            %    This will be the last 'n' logged variables to return. If
            %    the number of logged variables is less than this, all
            %    available data will be returned. Default is 1 if not
            %    supplied.
            %
            % Output
            %
            %  vals - The last 'n' logged values for the supplied variable,
            %    in the same shape as in the corresponding stored data
            %    field.
            %
            %
            % See Also: wsim.logger.logVal, wsim.logger.addVariable
            %
            
            if nargin < 3
                n = 1;
            end

            if ~isfield (obj.data, varname)
                error ('data logging field: %s does not exist', varname);
            end
            
            if obj.info.(varname).LastLogIndex < 1
                error ('No data has been logged for the variable %s yet.', varname);
            end
            
            % copy the pre-constructed indexing structure (created when
            % adding the variable)
            S = obj.info.(varname).IndexAssignment;
            
            % build the correct index into the logged variable by replacing
            % the appropriate index with the new log index
            if obj.info.(varname).LastLogIndex > n
                startindex = obj.info.(varname).LastLogIndex - n + 1;
            else
                startindex = 1;
            end
            
            S.subs{obj.info.(varname).IndexDimension} = startindex : obj.info.(varname).LastLogIndex;
            
            % get the values
            vals = subsref (obj.data.(varname), S);

        end
        
        
        function vals = getPrevLoggedVal(obj, varname, n)
            % get the (end - n) logged value of a variable
            %
            % Syntax
            %
            % val = getPrevLoggedVal(obj, varname, n)
            %
            % Description
            %
            % If a variable exists, the (last-n)th value logged to that
            % variable is returned.
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  varname - name of the variable to be for which to obtain the
            %    last logged value
            %
            %  n - the indicator of index of the logged variable to return.
            %   The value returned is the (last - n) logged value of the
            %   variable. i.e. to get the last value logged, use n = 0, the
            %   penultimate value n = 1 and so on.
            %
            % Output
            %
            %  vals - The last-n logged value for the supplied variable,
            %    in the same shape as in the corresponding stored data
            %    field.
            %
            %
            % See Also: wsim.logger.logVal, wsim.logger.addVariable
            %

            if ~isfield (obj.data, varname)
                error ('data logging field: %s does not exist', varname);
            end
            
            if obj.info.(varname).LastLogIndex < 1
                error ('No data has been logged for the variable %s yet.', varname);
            end
            
            % copy the pre-constructed indexing structure (created when
            % adding the variable)
            S = obj.info.(varname).IndexAssignment;
            
            % build the correct index into the logged variable by replacing
            % the appropriate index with the new log index
            logindex = obj.info.(varname).LastLogIndex - n;
            
            if logindex < 1
                error ('Requested value is before start of logging (or outside logging window)'); 
            end

            S.subs{obj.info.(varname).IndexDimension} = logindex;
            
            % get the value
            vals = subsref (obj.data.(varname), S);

        end
        
        
        function setSeries (obj, varname, newdata)
            % set values for an entire data series directly
            %
            % Syntax
            %
            % setSeries (obj, varname, newdata)
            %
            % Description
            %
            % setSeries set the values for an entire data series for a
            % logged variable directly. All existing logged data will be
            % wiped, and the new data will be put into the variable. If the
            % variable is associated with an independant variable, the
            % length of the new data series must match the length of the
            % associated independant variable.
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  varname - character vector with a single variable
            %   name for which the entire data series is to be set. If no
            %   exact match is found for the name in varname, a case
            %   insensitive search is performed and any unambiguous
            %   shortening of a variable name is allowed. e.g. to get the
            %   info for a logged variable named 'ExampleLoggedVariable'
            %   one can use the string 'ExampleLog'.
            %
            %  newdata - a new data series for the logged variable. This
            %   will replace any existing data.
            %
            %
            % See Also: wsim.logger.getInfo
            %

            assert (ischar (varname), 'varname should be a character vector containing the name of the variable for which to set the series.');
            
            % getInfo will return a structure with one field, which will be
            % the correct variable name (getInfo finds any unambiguous
            % match for the name)
            varinfo = obj.getInfo (varname);
            
            % there will only be one field in varinfo, which is the proper
            % name of the variable
            fname = fieldnames (varinfo);
            if ~strcmp(fname, varname)
                obj.logWarn (sprintf ('Using non-exact match in setSeries for %s. Input variable name was %s', fname, varname));
            end
            
            % replace the input variable name with the matched name
            varname = fname{1};
            
            % The field for this variable will contain the info for the
            % variable, so we get it's contents
            varinfo = varinfo.(varname);
            
            if ~isempty (varinfo.IndependentVariable)
                % check data length is same size as existing independent
                % variable
                ivinfo = obj.getInfo (varinfo.IndependentVariable);
                ivinfo = ivinfo.(varinfo.IndependentVariable);
                
                if size (newdata, varinfo.IndexDimension) ~= ivinfo.LastLogIndex ...
                    && size (newdata, varinfo.IndexDimension) ~= ivinfo.PreallocatedLogLength
                
                    error ('New data for series must be the same length as the linked independant variable');
                    
                end
                
            end
            
            % Check new data is the right size (matches existing variable
            % specification)
            sz = size (newdata);
            for ind = 1:numel (sz)
                if ind ~= varinfo.IndexDimension
                    if sz(ind) ~= varinfo.Size(ind)
                        error ('new data for variable series does not match existing variable size specification');
                    end
                end
            end
            
            % store the new data in the variable
            loggedvarsize = varinfo.LoggedSize;
            loggedvarsize(varinfo.IndexDimension) = varinfo.PreallocatedLogLength;
            % clear out existing data and replace with preallocated nans
            obj.data.(varname) = nan (loggedvarsize);
            
            % copy the pre-constructed indexing structure (created when
            % adding the variable)
            S = varinfo.IndexAssignment;
            
            % build the correct indices into the logged variable by
            % replacing the appropriate index with the new log indices
            S.subs{varinfo.IndexDimension} = 1:sz(varinfo.IndexDimension);
            
            % assign the new values to the series
            obj.data.(varname) = subsasgn (obj.data.(varname), S, newdata);

            % update the log index
            if ~isempty (varinfo.IndependentVariable)
                % assume the last correct log index is the same as the
                % linked independant variable
                obj.info.(varname).LastLogIndex = ivinfo.LastLogIndex;
            else
                % take last log index as end of provided series
                obj.info.(varname).LastLogIndex = sz(varinfo.IndexDimension);
            end

        end


        function setSilent(obj, bool)
            % turn off messages (but not warning or error messages)
            %
            % Syntax
            %
            % setSilent(obj, bool)
            %
            % Description
            %
            % Sets the silent property true or false. If true, then the
            % internal error messages of the logger class are not printed
            % out. Whatever the user specifies by setting the show variable
            % in the logWarn, logMesg functions are still printed out to
            % the command window.
            %
            %
            % Input
            %
            %  obj - wsim.logger object
            %
            %  bool - true/false flag. If true, messages will be silenced,
            %   but warning and error messages will still be displayed.
            %
            % See Also: 
            %
            
            check.isLogicalScalar (bool, true, 'bool');
            
            obj.silent = bool;
            
        end
        
        function info = getInfo (obj, reqvarnames)
            % get info on the logging settings for a variable
            %
            % Syntax
            %
            % info = getInfo (lg, reqvarnames)
            %
            % Description
            %
            % Returns a structure containing information on one or more
            % logged variables in the wsim.logger object. This is a subset
            % of the structure in the wsim.logger.info property. Variables
            % can be found by usng their full exact name or any unambiguous
            % cases insensitive matching shorter string.
            %
            % Input
            %
            %  lg - wsim.logger object
            %
            %  reqvarnames - character vector with a single variable
            %   name, or a cell array of character vectors with multiple
            %   variable names. If no exact match is found a case
            %   insensitive search is performed and any unambiguous
            %   shortening of a variable name is allowed. e.g. to get the
            %   info for a logged variable named 'ExampleLoggedVariable'
            %   one can use the string 'ExampleLog'.
            %
            % Output
            %
            %  info - structure array (one for each requested variable name
            %   in reqvarnames, with the following fields:
            %
            %    Description : string containing a description of the
            %     variable
            %
            %    Size : vector of integers specifying the size of the
            %      variable being logged for all dimensions.
            %
            %    PreallocatedLogLength : length of reserved log space for
            %      the variable
            %
            %    IndexAssignment : indexing structure used to aid indexing
            %      into the logged variable matrix (see help for subsasgn
            %      for more information)
            %
            %    IndexDimension : The matrix dimension along which the
            %      logging takes place
            %
            %    DataDimensions : A vector containing the dimesions which
            %      contain the variable data for each logged value of a
            %      variable. This excludes the log index dimension, and any
            %      singleton dimensions created in order to force logging
            %      along a certain dimension.
            %
            %    LastLogIndex : The index of the last logged item in the
            %      matrix
            %
            %    IndependentVariable : string containing the name of the
            %      independant variable associated with this variable. Can
            %      be empty if none is assigned, typically this is 'Time'.
            %
            %    LoggedSize : This is the size of the logged variable
            %      internally to the logger object, after adding any
            %      dimensions in order to log along the first singleton, or
            %      user specified dimension, but before preallocation.
            %
            %
            % See Also: wsim.logger.setSeries
            %

            if iscellstr (reqvarnames)
                % ok to proceed
            elseif ischar (reqvarnames)
                % make cell array
                reqvarnames = {reqvarnames};
            else
            	error ('field must be a char array, or cell array of strings');
            end
            
            varnames = fieldnames(obj.info);
            lvarnames = lower(varnames);
            
            for i = 1:numel(reqvarnames)
                
                if isfield (obj.info, reqvarnames{i})
                    % check for exact case-sensitive field match
                    info.(reqvarnames{i}) = obj.info.(reqvarnames{i});
                    
                else
                    % check for partial unambiguous match
                    
                    lreqvarname = lower(reqvarnames{i});

                    % check for case insensitive match of the start of the
                    % string
                    ind = find (strncmp (lreqvarname, lvarnames, length(lreqvarname)));

                    if isempty(ind)
                        error(['No matching logged variable name found for: ', reqvarnames{i}])
                    elseif length(ind) > 1
                        error('Ambiguous requested variable name: ''%s''. Could match one of %s', ...
                            reqvarnames{i}, sprintf ('%s ', varnames{ind}))
                    end
                    
                    % found one match, so copy it to the output variable
                    % and proceed
                    info.(varnames{ind}) = obj.info.(varnames{ind});

                end

            end
            
        end
        
        
        function truncateAllVariables (obj)
            % strips preallocated data space from all variables leaving only logged data
            %
            %
            
            fnames = fieldnames (obj.data);
            
            for ind = 1:numel (fnames)
                
                truncateVariable (obj, fnames{ind})
                
            end
            
        end
        
        
        function truncateVariable (obj, varname)
            % strips preallocated data space from one variable leaving only logged data
            %
            %
            
            assert (ischar (varname), 'varname must be a character vector');
            assert (isfield (obj.data, varname), ...
                'The copntents of varname does not appear to match exactly the name of any logged variables');
            
            if obj.info.(varname).PreallocatedLogLength > obj.info.(varname).LastLogIndex
                % copy the pre-constructed indexing structure (created when
                % adding the variable)
                S = obj.info.(varname).IndexAssignment;
            
                % build the correct index into the logged variable by replacing
                % the appropriate index with the new log index
                S.subs{obj.info.(varname).IndexDimension} = obj.info.(varname).LastLogIndex : obj.info.(varname).PreallocatedLogLength;
            
                % assign the new value
                obj.data.(varname) = subsasgn (obj.data.(varname), S, []);
            
                % update the preallocated length 
                obj.info.(varname).PreallocatedLogLength = obj.info.(varname).LastLogIndex;

            end
            
        end
        
    end % of methods
    
%     methods (Access=protected)
%         
%         function recPlotVar (obj, inds)
%             
%             for ind = 1:numel ()
%                 
%             end
%             
%         end
%         
%     end

end % of class.
