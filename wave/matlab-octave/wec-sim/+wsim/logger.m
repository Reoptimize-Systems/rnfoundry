classdef logger < handle
% wsim.logger
%
% Description
%
% A matlab class to create objects that can be used to log various events
% during the execution of a matlab script/function. This class is mainly a
% container class, with some additional functionality builtin to operate on
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
% (3) a global/persistent logger object can be created to log from various 
% functions
% (4) you can use matlab's event framework to log events by adding 
% appropriate listeners.
% 
% Examples
%
% Example 1
%
% 	l = wsim.logger;
%
% 	for i = 1:100,
% 		my_output_1 = 10*rand;
% 		height = 1.5*my_output_1 + 5*rand;
% 	
% 		l.logVal('weight', my_output_1);
% 		l.logIt(height);
% 	end
% 	
% 	l.setDesc('weight','Weight of Subjects');
% 	l.setDesc('height','Height of Subjects');
% 	l.setDefaultDesc('Subject ID');
% 	
% 	figure; l.plot2Vars('weight','height','LineWidth', 2, 'Color','r'); 
% 	figure; l.plotVar('weight','LineWidth', 2, 'Color','r'); 
% 	figure; l.plotVar('height','LineWidth', 2, 'Color','r'); 
% 	figure; l.plotFofVar('height',@log, 'LineWidth', 2, 'Color','r'); 
% 
% Also see logger_demo.m for example usage.
%
% Author: Richard Crozier
% Derived from the logger class developed by Pavan Mallapragada 
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
        errors; % Really? :)

        % Boolean value. If set to true, logger will not print any of its messages. It will still print the warnings/errors that the user specifies it to.
        silent;

        % The function handle used for plotting. Default is plot.
        plotfunc;

        % The function handle used for displaying interal messages of logger object. Default is set to warning.
        mesgfunc;

        % A string used for labeling the x-axis when plotting single variables.
        defaultDesc;
        
    end
    

    methods 
    % ==========================================================================
    % @brief constructs an empty logger object.
    % 
    % Initializes basic plotting functions, verbosity and messaging options.
    % Initializes the structures used for storing messages/warnings/errors
    % etc.
    %
    % @retval obj object of the logger class. 
    % ==========================================================================

        function obj = logger()
            % constructs an empty logger object.
            %
            % Syntax
            %
            % obj = wsim.logger ()
            %
            % Description
            %
            % Initializes basic plotting functions, verbosity and messaging
            % options. Initializes the structures used for storing
            % messages/warnings/errors etc.
            %
            % 

            obj.silent = false;

            obj.info = struct ();

            obj.messages = {};
            obj.warnings = {};
            obj.errors   = {};

            obj.plotfunc = @plot;
            obj.mesgfunc = @warning;
            obj.silent   = true;

            obj.defaultDesc = '';
            
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
            %    advance. 
            %
            %  'IndependentVariable' - string containing the name of
            %    another variable, already added to the logger, which is
            %    the independant variable for the new variable. Primarily
            %    used for plotting, but, also for preallocation. If the
            %    PreallocateStorage option is not specified, the new
            %    variable is by default preallocated the same storage as
            %    the independant variable.
            %

            
            options.Description = '';
            options.ForceLogDimension = [];
            options.PreallocateStorage = [];
            options.IndependentVariable = '';
            
            options = parse_pv_pairs (options, varargin);
            
            assert (isvarname (name), 'name must be a string containing a valid variable name');
            
            assert (isvector (varsize) && numel (varsize) >= 2, ...
                'varsize must be a vector of variable dimensions with two or more elements' );
            
            if ~isempty (options.PreallocateStorage)
                check.isPositiveScalarInteger (options.PreallocateStorage, true, 'PreallocateStorage');
            end
            
            assert (ischar (options.Description), 'Description must be a string');

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
                        indass.subs = {1, 1};
                        logdim = 1;
                    elseif numel (loggedvarsize) == 2 && loggedvarsize(2) == 1
                        indass.subs = {1, 1};
                        logdim = 2;
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
                                           'LastLogIndex', 0, ...
                                           'IndependentVariable', options.IndependentVariable, ...
                                           'LoggedVariableNumber', numel(obj.info) + 1);

            end
            
            % make the singleton logging dimension the required
            % length, the data will be preallocated with nans
            loggedvarsize(logdim) = options.PreallocateStorage;

            obj.data.(name) = nan (loggedvarsize);
            
            obj.fieldNames = fieldnames (obj.info);
            obj.numVariables = numel (obj.fieldNames);
            
        end

%     % ==========================================================================
%     % @brief generic logger function without a specific name for the logged-object.
%     % 
%     % This function accepts any matlab variable as input, and starts logging it
%     % with the same name.  This is useful when you do not want to give a different
%     % name to the logged object than its own variable name. This function
%     % determines if the variable is a numeric scalar or an object, and calls the
%     % appropriate named-logger function.
%     % 
%     % @param obj An object of the logger class.
%     % @param variable any variable from your workspace.
%     % ==========================================================================
% 
%         function [] = logData (obj, variable)
% 
%             if isnumeric(variable)
%                 obj.logVal(varname, variable);
%             else
%                 obj.logObj(varname, variable);
%             end
%             
%         end
 
    % ==========================================================================
    % @brief logs the given message in the error category. 
    %
    % Logs the error messages. This can be used for failed assertions or any other errors.
    % This is useful when a logger object is declared persistent or global and mutliple
    % functions are communicating with the same logger object.
    %
    % @param obj An object of the logger class.
    % @param mesg : a string containing your error message. If you are sure you will never set show to true, this can be any object. Although it is not recommended.
    % @param show : true/false specifying whether to print the message on the command window.
    % ==========================================================================

        function [] = logErr(obj, mesg, show)
            obj.errors{end+1} = mesg;
            if nargin > 2 && show
                fprintf('(e): %s\n',mesg);
            end
        end
 
    % ==========================================================================
    % @brief logs the given message in the warning category. 
    %
    % Logs the warning messages. 
    %
    % @param obj An object of the logger class.
    % @param mesg : a string containing your warning message. If you are sure you will never set show to true, this can be any object. Although it is not recommended.
    % @param show : true/false specifying whether to print the message on the command window.
    % ==========================================================================
        function [] = logWarn(obj, mesg, show)
            obj.warnings{end+1} = mesg;
            if nargin > 2 && show
                fprintf('(w): %s\n',mesg);
            end
        end
 
    % ==========================================================================
    % @brief logs the given message in the information category. 
    %
    % Logs any information messages given. 
    %
    % @param obj An object of the logger class.
    % @param mesg : a string containing your warning message. If you are sure you will never set show to true, this can be any object. Although it is not recommended.
    % @param show : true/false specifying whether to print the message on the command window.
    % ==========================================================================

        function [] = logMesg(obj, mesg, show)
            obj.messages{end+1} = mesg;
            if nargin > 2 && show
                fprintf('(i): %s\n',mesg);
            end
        end

    % ==========================================================================
    % @brief This function logs the non-numeric objects, and stores them in a field as a cell array.
    %
    % This is useful for storing arrays, structures, etc that are non-scalar and non-numeric. These fields cannot be used for
    % plotting at this point. More functionality on this would be added in the newer versions. At this point
    % The class is only a container for these objects. 
    % 
    % @param obj An object of the logger class.
    % @param field The field name you are logging.
    % @param val   The value any matlab object.
    % ==========================================================================
        function [] = logObj(obj, field, val)
            try
                obj.data.(field){end+1} = val;
            catch
                obj.addprop(field);
                if ~obj.silent
                    fprintf(['Logger: new field added to the logger object: ', field,'\n']);
                end
                obj.fieldNames{end+1} = field;
                obj.data.(field) = {};
                obj.data.(field){end+1} = val;
            end

        end

    % ==========================================================================
    % @brief This function sets the default message function for the logger object.
    %
    % If you want logger to not quit your process, because of an error it
    % encounters, you can set it to @warning (it is default). If you want to be
    % strict, and debut errors in logger function calls, you can set it to @error.
    % 
    % @param obj An object of the logger class.
    % @param mf function handle you would like to use for displaying logger's internal error/warning messages.
    % ==========================================================================
        function [] = setMesgFunc(obj, mf) 
            obj.mesgfunc  = mf; 
        end 

    % ==========================================================================
    % @brief Default description which is used for labeling the x-axis.
    % 
    % This is useful while plotting single-variables (fields).
    % ==========================================================================
        function [] = setDefaultDesc(obj, str)
            obj.defaultDesc = str;
        end
 
    % ==========================================================================
    % @brief This function is not used yet. This is written here for possible future use.
    % ==========================================================================
        function [] = setDefaultField(obj,field)
            obj.defaultField =  field;
        end

    % ==========================================================================
    % @brief Overriden disp function for the logger class.
    %
    % Displays the logger object with all its fields and their sizes.
    % ==========================================================================

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

    % ==========================================================================
    % @brief set the description for a field to be used for x-y labels in plotting. 
    % 
    % The arguments could be two cell arrays, with fieldname
    % and descritpion corresponding to each other, or two strings. 
    %
    % @param obj An object of the logger class
    % @param f   A string or a cell array of strings containing the field names 
    % @param desc   A string or a cell array of strings containing descriptions of corresponding fields
    %
    % ==========================================================================

        function [obj] = setDesc(obj, field, desc)

            if ischar(field) && ischar(desc)
                obj.info.(field).Description = desc;
            elseif iscell(field) && iscell(desc)
                
                assert (samesize (field, desc), ...
                    'If field and desc are cell string arrays they must be the same size');
                
                for i = 1:numel(field)
                    obj.info.(field{i}).Description = desc{i};
                end
            else
                error ('field and desc must be both either strings, or cell string arrays of the same size')
            end
            
        end

    % ==========================================================================
    % @brief Sets the plotting function to be used. 
    % 
    % Sets the plot function. Its default vaule is plot. You can make it semilogx, semilogy, etc.
    %
    % @param pf_handle Plot function handle.
    % ==========================================================================

        function [obj] = setPlotFunc(obj, pf_handle)
            obj.plotfunc = pf_handle;
        end

    % ==========================================================================
    % @brief given a string specifying a numeric scalar field, and a funtion handle,
    % the function is first applied to the field, and then it is plotted.
    %  
    % For example, if a user logs a fields 'height', the log(heights) can be plotted as
    % logger.plotFofVar('height',@log). The parameters to
    % the plot can be provided after the fields, and all
    % those arguments go to the plot function.  For example,
    % logger.plotFofVar('height','@log','LineWidth',2,'Color','r'); will pass the
    % last four arguments to plot.
    %
    % @param obj logger object
    % @param field   a string specifying the name of logged field.
    % @param func   a function handle that is to be evaluated on the field.
    % @param varargin  any arguments that are to be sent to the plotting function.
    % 
    % @retval h A handle to the plot generated. Useful for formatting by the user.
    % ==========================================================================


        function [h] = plotFofVar(obj, field, func, varargin)

            if ~ischar(field), 	 mesgfunc([field 'must be a string specifying field that are already added to the logger object.']); return; end
            if ~isnumeric(obj.data.(field)) mesgfunc(['Plotting only numeric values is supported at this point. Not generating the plot for' field]); return; end
            if ~isa(func, 'function_handle') mesgfunc(['third argument func must be a function handle']); return; end

            plotvals = func(obj.data.(field));

            h = obj.plotfunc(obj.data.(field), varargin{:});

            hold on;

            ylabel([func2str(func) ' (' obj.info.Descriptions.(field) ')'], 'FontSize', 16);
            if ~isempty(obj.defaultDesc)
                xlabel(obj.defaultDesc,'FontSize',16);
            end

            set(gca,'FontSize',16);

            hold off;
        end

    % ==========================================================================
    % @brief given a string specifying a numeric scalar field, it is plotted.
    %  
    % For example, if a user logs a fields 'height', using
    % logger.plotVar('height') will plot height. The parameters to
    % the plot can be provided after the fields, and all
    % those arguments go to the plot function.  For example,
    % logger.plotVar('height','LineWidth',2,'Color','r'); will pass the
    % last four arguments to plot.
    %
    % @param obj logger object
    % @param field   a string specifying the name of logged field 1.
    % @param varargin  any arguments that are to be sent to the plotting function.
    % 
    % @retval h A handle to the plot generated. Useful for formatting by the user.
    % ==========================================================================


        function [h] = plotVar(obj, field, varargin)

            if ~ischar(field), 	 mesgfunc([field 'must be a string specifying field that are already added to the logger object.']); return; end
            if ~isnumeric(obj.data.(field)) mesgfunc(['Plotting only numeric values is supported at this point. Not generating the plot for' field]); return; end


            h = obj.plotfunc(obj.data.(field), varargin{:});

            hold on;

            ylabel(obj.info.Descriptions.(field), 'FontSize', 16);

            if ~isempty(obj.defaultDesc)
                xlabel(obj.defaultDesc,'FontSize',16);
            end

            set(gca,'FontSize',16);
        end


        % ==========================================================================
        % @brief given two numeric scalar fields, they are plotted one against another. 
        %  
        %  
        % For example, if a user logs two fields 'height' and 'weight', giving
        % logger.plotvars('height','weight') will plot height vs. weight. The parameters to
        % the plot can be provided after the fields, and all
        % those arguments go to the plot function.  For example,
        % logger.plotvars('height','weight','LineWidth',2,'Color','r'); will pass the
        % last four arguments to plot.
        %
        % @param obj logger object
        % @param f1   a string specifying the name of logged field 1.
        % @param f2   a string specifying the name of logged field 1.
        % @param varargin  any arguments that are to be sent to the plotting function.
        % 
        % @retval h A handle to the plot generated. Useful for formatting by the user.
        % ==========================================================================

        function [h] = plot2Vars(obj, f1, f2, varargin)
            
            options.PlotFcnArgs = {};
            
            options = parse_pv_pairs (options, varargin);

            if ~ischar(f1) || ~ischar(f2)
                feval (obj.mesgfunc, 'f1 and f2 must be strings specifying fields that are already added to the logger object.');
            end

            if ~isnumeric(obj.(f1)) || ~isnumeric(obj.(f2))
                feval (obj.mesgfunc, sprintf ('Plotting only numeric values is supported at this point. Not generating the plot %s vs %s', f1, f2));
            end
            
            if ( numel (obj.info.(f1).Size) ~= numel (obj.info.(f2).Size) ) ...
                feval (obj.mesgfunc, 'Cannot plot variables of different sizes against each other');
            end

            legstrs = {};
            
            h = obj.plotfunc(obj.(f1), obj.(f2), options.PlotFcnArgs{:});

            hold on;

            desc1 = obj.info.Descriptions.(f1);
            desc2 = obj.info.Descriptions.(f2);

            xlabel(desc1, 'FontSize', 16);
            ylabel(desc2, 'FontSize', 16);
            set(gca,'FontSize',16);
            
        end


        function status = logVal(obj, varname, val, ignoremissing)
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
            % See Also: 
            %

            status = 0;
            
            if nargin < 4
                ignoremissing = false;
            end
            
            if ~isfield (obj.data, varname)
                if ignoremissing
                    status = -1;
                    return;
                else
                    error ('data logging field: %s does not exist', varname);
                end
            end
            
            % copy the pre-constructed indexing structure (created when
            % adding the variable)
            S = obj.info.(varname).IndexAssignment;
            
            % build the correct index into the logged variable by replacing
            % the appropriate index with the new log index
            S.subs{obj.info.(varname).IndexDimension} = obj.info.(varname).LastLogIndex + 1;
            
            % assign the new value
            obj.data.(varname) = subsasgn (obj.data.(varname), S, val);
            
            % increment the data index counter for this field
            obj.info.(varname).LastLogIndex = obj.info.(varname).LastLogIndex + 1;

        end
        
        function setSeries (obj, varname, newdata)
            % set values for an entire data series directly (replaces
            % existing data)
            
            varinfo = obj.getInfo (varname);
            
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
            obj.data.(varname) = newdata;
            
        end


        function setSilent(obj, bool)
            % turn off messages (but not warining or errors messages)
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
            % Output
            %
            %
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
            % logged variable in the wsim.logger object. This is a subset
            % of the structure in the wsim.logger.info property. Variables
            % can be found by usng their full exact name or any unambiguous
            % cases insensitive matching shorter string.
            %
            % Input
            %
            %  lg - wsim.logger object
            %
            %  reqvarnames - string (char array) with a single variable
            %    name, or a cell array of strings with multiple variable
            %    names. If no exact match is found a case insensitive
            %    search is performed and any unambiguous shortening of a
            %    variable name is allowed. e.g. to get the info for a
            %    logged variable named 'ExampleLoggedVariable' one can use
            %    the string 'ExampleLog'.
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
            %    LastLogIndex : The index of the last logged item in the
            %      matrix
            %
            %    IndependentVariable : string containing the name of the
            %      independant variable associated with this variable. Can
            %      be empty if none is assigned, typically this is 'Time'.
            %
            %
            % See Also: 
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
        
    end % of methods

end % of class.
