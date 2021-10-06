function setaddaxisprops(axnum, varargin)
% setaddaxisprops: sets the properties of a given axis added by the addaxis
% function
%
% Syntax
%
% setaddaxisprops(axnum, 'PropertyName', 'PropertyValue', ...)
%
% Input
%
%   axnum - axis number you wish to change the properties of or, a string
%           'all'. If 'all' is passed in the proerties are applied to all
%           axes.
%
%   The pv-paris following axnum are the same that can be specisfied in
%   'axes'.
    
    %  store the current/main axes handle
    cah = gca;

    % get the handles of all the axes the first of these ought to be the
    % main axes
    axh = getaddaxisdata(cah,'axisdata');

    % Loop through all the main axes getting the subaxes for each
    axhand = cah;
    for i = 1:length(axh)
        axhand(i+1) = axh{i}(1);
    end

    if ischar(axnum)

        if strcmp(axnum,'all')

            for i = 1:length(axhand)
                
                %  set current axis to the axis to be modified
                axes(axhand(i));

                % Set the properties passed in as pv-pairs
                set(gca, varargin{:});
                
            end
            
            %  set current axis back to the original value
            axes(cah);
        else
            %  set current axis back to the original value
            axes(cah);
            % throw error
            error('Invalid string entered, must be ''all'' or an integer specifying the axis number')
        end

    else
        
        %  set current axis to the axis to be modified
        axes(axhand(axnum));

        % Set the properties passed in as pv-pairs
        set(gca, varargin{:});

        %  set current axis back to the original value
        axes(cah);

    end

end