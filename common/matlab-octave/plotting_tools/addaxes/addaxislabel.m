function varargout = addaxislabel(axnum, varargin)
%ADDAXISLABEL adds axis labels to axes made with ADDAXIS.m
%
% Syntax
%
%  handle_to_text = addaxislabel(axis_number, label, 'PropertyName', PropertyValue', ...);
%
%  See also
%  ADDAXISPLOT, ADDAXIS, SPLOT 
  
    % get current axes handle
    cah = gca;
    
    % get the handles to all the axes in the figure
    axh = getaddaxisdata(cah,'axisdata');

    %  get axis handles to all the main axes and subaxes in the order they
    %  were added
    axhand = cah;
    postot(1,:) = get(cah,'position');
    for I = 1:length(axh)
        axhand(I+1) = axh{I}(1);
        postot(I+1,:) = get(axhand(I+1),'position');
    end

    %  set current axis to the axis to be labeled
    axes(axhand(axnum));
    
    % Set the label of the y axis and get the handle to the label object so
    % we can change it's colour
    htxt = ylabel(varargin{:});
    
    % Set the label text colour to the same as the axis colour
    set(htxt,'color',get(axhand(axnum),'ycolor'));

    %  set current axis back to the main axis
    axes(cah);

    if nargout == 1
        % return the handle to the label object if requested
        varargout{1} = htxt;
    end

end
