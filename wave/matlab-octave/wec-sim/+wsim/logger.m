% ==========================================================================
%> @file logger.m
%> @brief MATLAB Logger
% ==========================================================================
% ==========================================================================
%
%> A matlab class to create objects that can be used to log various events during the
%> exectuion of a matlab script/function. This class is mainly a container class, with
%> some additional functionality builtin to operate on the stored data in generating
%> plots, applying functions etc.
% 
%> The main aim of this class is to consolidate and store all the outputs and
%> messages from a MATLAB script/function into a single object.  These objects
%> can be saved, and retrieved later, with all the information in one place
%> Several utility functions (mainly for plotting) make it easy to operate on
%> the stored data.
%
%> There are several ways in which you can use this class. You can either
%> (1) create a logger object, and start logging into the class
%> (2) user's class can be inherited from logger
%> (3) a global/persistent logger object can be created to log from various functions
%> (4) you can use matlab's event framework to log events by adding appropriate listeners.
% 
%> Simple Example usage:
%> 
%> 	l = logger;
%>
%> 	for i = 1:100,
%> 		my_output_1 = 10*rand;
%> 		height = 1.5*my_output_1 + 5*rand;
%> 	
%> 		l.logVal('weight', my_output_1);
%> 		l.logIt(height);
%> 	end
%> 	
%> 	l.setDesc('weight','Weight of Subjects');
%> 	l.setDesc('height','Height of Subjects');
%> 	l.setDefaultDesc('Subject ID');
%> 	
%> 	figure; l.plot2Vars('weight','height','LineWidth', 2, 'Color','r'); 
%> 	figure; l.plotVar('weight','LineWidth', 2, 'Color','r'); 
%> 	figure; l.plotVar('height','LineWidth', 2, 'Color','r'); 
%> 	figure; l.plotFofVar('height',@log, 'LineWidth', 2, 'Color','r'); 
%> 
%> Also see logger_demo.m for example usage.
%
% ==========================================================================
% Author: 	 Pavan Mallapragada 
% Organization:  Massachusetts Institute of Technology
% Contact:       <pavan_m@mit.edu>
% Created:       Oct 01, 2011 
% ==========================================================================

classdef logger < handle

    properties (GetAccess = public, SetAccess = private)
        
        %> A structure that stores the descriptions of the fields. These descriptions are used for labeling plots.
        info;
        
        %> Number of variables being logged
        numVariables;

        %> A cell array that stores all the fields that are being logged by the object currently, as strings.
        fieldNames;

        %> A structure in which the logged data will be stored.
        data;

        %> A string containing field name to perform default operations on.  Currently this is unused.
        defaultfield;

        %> Cell array to store messages
        messages;

        %> Cell array to store warnings
        warnings;

        %> Cell array to store error messages
        errors; % Really? :)

        %> Boolean value. If set to true, logger will not print any of its messages. It will still print the warnings/errors that the user specifies it to.
        silent;

        %> The function handle used for plotting. Default is plot.
        plotfunc;

        %> The function handle used for displaying interal messages of logger object. Default is set to warning.
        mesgfunc;

        %> A string used for labeling the x-axis when plotting single variables.
        defaultDesc;
        
    end
    

    methods 
    % ==========================================================================
    %> @brief constructs an empty logger object.
    %> 
    %> Initializes basic plotting functions, verbosity and messaging options.
    %> Initializes the structures used for storing messages/warnings/errors
    %> etc.
    %>
    %> @retval obj object of the logger class. 
    % ==========================================================================

        function [obj] = logger()

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
            
            options.Description = '';
            options.ForceLogDimension = [];
            options.PreallocateStorage = 1;
            options.IndependentVariable = '';
            
            options = parse_pv_pairs (options, varargin);
            
            assert (isvarname (name), 'name must be a string containing a valid variable name');
            
            assert (isvector (varsize) && numel (varsize) >= 2, ...
                'varsize must be a vector of variable dimensions with two or more elements' );
            
            check.isPositiveScalarInteger (options.PreallocateStorage, true, 'PreallocateStorage');
            
            assert (ischar (options.Description), 'Description must be a string');

            if ~isempty (options.IndependentVariable)
                
                assert (isvarname (options.IndependentVariable), ...
                    'IndependentVariable must be a string containing a valid variable name');
                
                assert (isfield (obj.info, options.IndependentVariable), [...
'If specifiying an independent variable for a logged variable, you must add\n' ...
'the independent variable first. The specified independent variable, "%s", \n' ...
'was not found in this logger object.' ]);

            end
            
            % generate input structure for subasgn function, adding
            % to log is done using this function to add data along the
            % first dimesion with length 1
            indass.type = '()';
            logdim = 0;
            
            if ~isempty (options.ForceLogDimension)
                
                check.isPositiveScalarInteger (options.ForceLogDimension, true, 'ForceLogDimension');
                
                % use user specified dimension
                if options.ForceLogDimension > numel (varsize)
                    % add the necessary singleton dimensions to varsize
                    for ind = (numel (varsize) + 1):options.ForceLogDimension
                        varsize(ind) = 1;
                    end
                else
                    % check this dimension is of size 1
                    assert (varsize(options.ForceLogDimension) == 1, ...
                        'The size of data dimension chosen in ForceLogDimension must be 1');
                end
                
                indass.subs = repmat ({':'}, 1, numel (varsize));
                indass.subs{options.ForceLogDimension} = 1;
                logdim = options.ForceLogDimension;
                    
            else
                % find the first singleton dimension of array and log along
                % it
                
                for ind = 1:numel (varsize)
                    if varsize(ind) == 1
                        logdim = ind;
                    end
                end
                
                if logdim == 0
                    % no singleton dimensions were found, so add one at the
                    % end
                    varsize(end+1) = 1;
                    logdim = numel (varsize);
                    indass.subs = repmat ({':'}, 1, numel (varsize));
                else
                    if numel (varsize) == 2 && varsize(1) == 1
                        indass.subs = {1, 1};
                        logdim = 1;
                    elseif numel (varsize) == 2 && varsize(2) == 1
                        indass.subs = {1, 1};
                        logdim = 2;
                    else
                        indass.subs = repmat ({':'}, 1, numel (varsize));
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
                                           'IndependentVariable', options.IndependentVariable );

            end
            
            % make the singleton logging dimension the required
            % length, the data will be preallocated with nans
            varsize(logdim) = options.PreallocateStorage;

            obj.data.(name) = nan (varsize);
            
            obj.fieldNames = fieldnames (obj.info);
            obj.numVariables = numel (obj.fieldNames);
            
        end

%     % ==========================================================================
%     %> @brief generic logger function without a specific name for the logged-object.
%     %> 
%     %> This function accepts any matlab variable as input, and starts logging it
%     %> with the same name.  This is useful when you do not want to give a different
%     %> name to the logged object than its own variable name. This function
%     %> determines if the variable is a numeric scalar or an object, and calls the
%     %> appropriate named-logger function.
%     %> 
%     %> @param obj An object of the logger class.
%     %> @param variable any variable from your workspace.
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
    %> @brief logs the given message in the error category. 
    %>
    %> Logs the error messages. This can be used for failed assertions or any other errors.
    %> This is useful when a logger object is declared persistent or global and mutliple
    %> functions are communicating with the same logger object.
    %>
    %> @param obj An object of the logger class.
    %> @param mesg : a string containing your error message. If you are sure you will never set show to true, this can be any object. Although it is not recommended.
    %> @param show : true/false specifying whether to print the message on the command window.
    % ==========================================================================

        function [] = logErr(obj, mesg, show)
            obj.errors{end+1} = mesg;
            if nargin > 2 && show
                fprintf('(e): %s\n',mesg);
            end
        end
 
    % ==========================================================================
    %> @brief logs the given message in the warning category. 
    %>
    %> Logs the warning messages. 
    %>
    %> @param obj An object of the logger class.
    %> @param mesg : a string containing your warning message. If you are sure you will never set show to true, this can be any object. Although it is not recommended.
    %> @param show : true/false specifying whether to print the message on the command window.
    % ==========================================================================
        function [] = logWarn(obj, mesg, show)
            obj.warnings{end+1} = mesg;
            if nargin > 2 && show
                fprintf('(w): %s\n',mesg);
            end
        end
 
    % ==========================================================================
    %> @brief logs the given message in the information category. 
    %>
    %> Logs any information messages given. 
    %>
    %> @param obj An object of the logger class.
    %> @param mesg : a string containing your warning message. If you are sure you will never set show to true, this can be any object. Although it is not recommended.
    %> @param show : true/false specifying whether to print the message on the command window.
    % ==========================================================================

        function [] = logMesg(obj, mesg, show)
            obj.messages{end+1} = mesg;
            if nargin > 2 && show
                fprintf('(i): %s\n',mesg);
            end
        end

    % ==========================================================================
    %> @brief This function logs the non-numeric objects, and stores them in a field as a cell array.
    %>
    %> This is useful for storing arrays, structures, etc that are non-scalar and non-numeric. These fields cannot be used for
    %> plotting at this point. More functionality on this would be added in the newer versions. At this point
    %> The class is only a container for these objects. 
    %> 
    %> @param obj An object of the logger class.
    %> @param field The field name you are logging.
    %> @param val   The value any matlab object.
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
    %> @brief This function sets the default message function for the logger object.
    %>
    %> If you want logger to not quit your process, because of an error it
    %> encounters, you can set it to @warning (it is default). If you want to be
    %> strict, and debut errors in logger function calls, you can set it to @error.
    %> 
    %> @param obj An object of the logger class.
    %> @param mf function handle you would like to use for displaying logger's internal error/warning messages.
    % ==========================================================================
        function [] = setMesgFunc(obj, mf) 
            obj.mesgfunc  = mf; 
        end 

    % ==========================================================================
    %> @brief Default description which is used for labeling the x-axis.
    %> 
    %> This is useful while plotting single-variables (fields).
    % ==========================================================================
        function [] = setDefaultDesc(obj, str)
            obj.defaultDesc = str;
        end
 
    % ==========================================================================
    %> @brief This function is not used yet. This is written here for possible future use.
    % ==========================================================================
        function [] = setDefaultField(obj,field)
            obj.defaultField =  field;
        end

    % ==========================================================================
    %> @brief Overriden disp function for the logger class.
    %>
    %> Displays the logger object with all its fields and their sizes.
    % ==========================================================================

        function [] = disp(obj)
            fprintf('Logger Object\n\n');
            fprintf('User Fields:\n');

            for i = 1:numel(obj.fieldNames)
                fprintf('%12s : %d log entries\n', obj.fieldNames{i}, length(obj.(obj.fieldNames{i})));
            end

            fprintf('OTHERS:\n');

            fprintf('    messages : %d log entries\n', length(obj.messages));
            fprintf('    warnings : %d log entries\n', length(obj.warnings));
            fprintf('      errors : %d log entries\n', length(obj.errors));
        end

    % ==========================================================================
    %> @brief set the description for a field to be used for x-y labels in plotting. 
    %> 
    %> The arguments could be two cell arrays, with fieldname
    %> and descritpion corresponding to each other, or two strings. 
    %>
    %> @param obj An object of the logger class
    %> @param f   A string or a cell array of strings containing the field names 
    %> @param desc   A string or a cell array of strings containing descriptions of corresponding fields
    %>
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
    %> @brief Sets the plotting function to be used. 
    %> 
    %> Sets the plot function. Its default vaule is plot. You can make it semilogx, semilogy, etc.
    %>
    %> @param pf_handle Plot function handle.
    % ==========================================================================

        function [obj] = setPlotFunc(obj, pf_handle)
            obj.plotfunc = pf_handle;
        end

    % ==========================================================================
    %> @brief given a string specifying a numeric scalar field, and a funtion handle,
    %> the function is first applied to the field, and then it is plotted.
    %>  
    %> For example, if a user logs a fields 'height', the log(heights) can be plotted as
    %> logger.plotFofVar('height',@log). The parameters to
    %> the plot can be provided after the fields, and all
    %> those arguments go to the plot function.  For example,
    %> logger.plotFofVar('height','@log','LineWidth',2,'Color','r'); will pass the
    %> last four arguments to plot.
    %>
    %> @param obj logger object
    %> @param field   a string specifying the name of logged field.
    %> @param func   a function handle that is to be evaluated on the field.
    %> @param varargin  any arguments that are to be sent to the plotting function.
    %> 
    %> @retval h A handle to the plot generated. Useful for formatting by the user.
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
    %> @brief given a string specifying a numeric scalar field, it is plotted.
    %>  
    %> For example, if a user logs a fields 'height', using
    %> logger.plotVar('height') will plot height. The parameters to
    %> the plot can be provided after the fields, and all
    %> those arguments go to the plot function.  For example,
    %> logger.plotVar('height','LineWidth',2,'Color','r'); will pass the
    %> last four arguments to plot.
    %>
    %> @param obj logger object
    %> @param field   a string specifying the name of logged field 1.
    %> @param varargin  any arguments that are to be sent to the plotting function.
    %> 
    %> @retval h A handle to the plot generated. Useful for formatting by the user.
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
        %> @brief given two numeric scalar fields, they are plotted one against another. 
        %>  
        %>  
        %> For example, if a user logs two fields 'height' and 'weight', giving
        %> logger.plotvars('height','weight') will plot height vs. weight. The parameters to
        %> the plot can be provided after the fields, and all
        %> those arguments go to the plot function.  For example,
        %> logger.plotvars('height','weight','LineWidth',2,'Color','r'); will pass the
        %> last four arguments to plot.
        %>
        %> @param obj logger object
        %> @param f1   a string specifying the name of logged field 1.
        %> @param f2   a string specifying the name of logged field 1.
        %> @param varargin  any arguments that are to be sent to the plotting function.
        %> 
        %> @retval h A handle to the plot generated. Useful for formatting by the user.
        % ==========================================================================

        function [h] = plot2Vars(obj, f1, f2, varargin)

            if ~ischar(f1) || ~ischar(f2) 
                mesgfunc([f1 ' and ' f2 'must be strings specifying fields that are already added to the logger object.']);
            end


            if ~isnumeric(obj.(f1)) || ~isnumeric(obj.(f2))
                mesgfunc(['Plotting only numeric values is supported at this point. Not generating the plot' f1 'vs' f2]);
            end

            h = obj.plotfunc(obj.(f1), obj.(f2),varargin{:});

            hold on;

            desc1 = obj.info.Descriptions.(f1);
            desc2 = obj.info.Descriptions.(f2);

            xlabel(desc1, 'FontSize', 16);
            ylabel(desc2, 'FontSize', 16);
            set(gca,'FontSize',16);
        end


        
        % ==========================================================================
        %> @brief A function that logs a value in a Matlab array, specified by its name.
        %>
        %> If a field exists, the element is added to the next index of the 
        %> field.
        %> 
        %> @param obj An object of the logger class.
        %> @param field the name of the field to store the logged values in.
        %> @param value : numeric value to be logged. Must be of the size
        %>        specified when the variable was added (usinig addVariable)
        %> @param ignoremissing flag determining whether to throw an error
        %>        if field is not found, or just silently return without taking any
        %>        action
        %>
        %> @retval status A flag indicating if the operation was successful
        % ==========================================================================

        function status = logVal(obj, field, val, ignoremissing)
            
            status = 0;
            
            if nargin < 4
                ignoremissing = false;
            end
            
            if ~isfield (obj.data, field)
                if ignoremissing
                    status = -1;
                    return;
                else
                    error ('data logging field: %s does not exist', field);
                end
            end
            
            % copy the pre-constructed indexing structure (created when
            % adding the variable)
            S = obj.info.(field).IndexAssignment;
            
            % build the correct index into the logged variable by replacing
            % the appropriate index with the new log index
            S.subs{obj.info.(field).IndexDimension} = obj.info.(field).LastLogIndex + 1;
            
            % assign the new value
            obj.data.(field) = subsasgn (obj.data.(field), S, val);
            
            % increment the data index counter for this field
            obj.info.(field).LastLogIndex = obj.info.(field).LastLogIndex + 1;

        end

        % ==========================================================================
        %> @brief To set true or false for 'silent' behavior; No messages from the logger class.
        %>
        %> This function sets the silent variable true or false. If true, then the internal error messages
        %> of the logger class are not printed out. Whatever the user specifies by setting the show variable
        %> in the logWarn, logMesg functions are still printed out to the command window.
        %>
        %> @param obj An object of the logger class.
        %> bool : true or false specifying the silent behavior of the class.
        % ==========================================================================
        function [] = setSilent(obj, bool)
            
            check.isLogicalScalar (bool, true, 'bool');
            
            obj.silent = bool;
            
        end
        
    end % of methods

end % of class.
