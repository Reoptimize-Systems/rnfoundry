function varargout = addaxislabel(input, varargin)
%ADDAXISLABEL adds axis labels to axes made with ADDAXIS.m
%
% Syntax
%
%  handle_to_text = addaxislabel(axis_number, label, 'PropertyName', PropertyValue', ...);
%
%  See also
%  ADDAXISPLOT, ADDAXIS, SPLOT 

    assert (isscalar (input), 'input must be a scalar integer, or single axes handle');
  
    ret = what_is_input (input);
    
    switch ret
        
        case 'label'
            
            % get current axes handle
            current_ax = gca;
            axnum = input;
            
            % get the handles to all the axes in the figure
            axh = getaddaxisdata (current_ax, 'axisdata');

            %  get axis handles to all the main axes and subaxes in the order they
            %  were added
            axhand = current_ax;

            postot(1,:) = get (current_ax, 'position');

            for I = 1:length(axh)
                axhand(I+1) = axh{I}(1);
                postot(I+1,:) = get (axhand(I+1), 'position');
            end
            
            change_ax = axhand(axnum);

        case 'axes_handle'
            
            change_ax = input;
            
        otherwise
            
            error ('Invalid input, input must be the number of the axes, or axes handle for the axes for which the ylabel is to be set.');
            
    end
    
    % Set the label of the y axis and get the handle to the label object so
    % we can change it's colour
    htxt = ylabel (change_ax, varargin{:});
    
    % Set the label text colour to the same as the axis colour
    set ( htxt, 'color', get (change_ax, 'ycolor') );

    if nargout == 1
        % return the handle to the label object if requested
        varargout{1} = htxt;
    end

end

function ret = what_is_input (thing)

    if check.isAxes (thing, false)
        ret = 'axes_handle';
    elseif check.isint2eps (thing)
        ret = 'label';
    else
        ret = '';
    end

end

