function output_txt = addaxisdtips(obj, event_obj)
% Display an observation's Y-data and label for a datatip
% obj          Currently not used (empty)
% event_obj    Handle to event object
%
% output_txt   Datatip text (string or string cell array)
%
% This datacursor callback calculates the actual data in a plot which has
% been created using addaxis, and hence has it's y data rescaled. 
%

    pos = get(event_obj, 'Position');
    x = pos(1); y = pos(2);

    % The information necessary to reconstruct the actual data from the
    % position is stored in a field named 'plotdata' in the plot's
    % application data. This is set by aa_splot.
    plotdata = getappdata(get(event_obj, 'Target'), 'plotdata');

    % Do some checking then work out what the original data was using the
    % supplied y axis values
    if isstruct(plotdata) && all(isfield(plotdata, {'yl2', 'yl'}))
        y = ((y - plotdata.yl(1)) ./ (plotdata.yl(2) - plotdata.yl(1))) ...
            .* (plotdata.yl2(2) - plotdata.yl2(1)) + plotdata.yl2(1);
    end

    % Set the output text for the datatip
    output_txt = {['X: ', num2str(x,4)], ['Y: ',num2str(y,4)]};

end



